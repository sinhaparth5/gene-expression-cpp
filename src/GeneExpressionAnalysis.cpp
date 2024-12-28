#include "GeneExpressionAnalysis.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iomanip>

namespace geneexp {

GeneExpressionAnalysis::GeneExpressionAnalysis()
    : isDataLoaded(false), isPreprocessed(false) {
    state.isNormalized = false;
    state.isLogTransformed = false;
    state.isFiltered = false;
}

GeneExpressionAnalysis::~GeneExpressionAnalysis() = default;

bool GeneExpressionAnalysis::loadData(const std::string& expressionFile, 
                                    const std::string& labelsFile) {
    try {
        // Load expression data
        std::ifstream expFile(expressionFile);
        if (!expFile.is_open()) {
            throw DataLoadingError("Cannot open expression file: " + expressionFile);
        }

        // Read header
        std::string line;
        std::getline(expFile, line);
        std::istringstream iss(line);
        std::string token;
        
        // Parse header to get sample labels
        while (std::getline(iss, token, ',')) {
            if (token != "gene_id") // Skip first column header
                data.sampleLabels.push_back(token);
        }

        // Read expression data
        std::vector<std::vector<double>> tempData;
        while (std::getline(expFile, line)) {
            std::vector<double> row;
            std::istringstream iss(line);
            
            // Get gene name
            std::getline(iss, token, ',');
            data.geneNames.push_back(token);
            
            // Get expression values
            while (std::getline(iss, token, ',')) {
                try {
                    row.push_back(std::stod(token));
                } catch (const std::exception& e) {
                    addWarning("Invalid expression value found: " + token);
                    row.push_back(0.0); // Default to zero for invalid values
                }
            }
            tempData.push_back(row);
        }

        // Convert to Eigen matrix
        data.expressionMatrix = Eigen::MatrixXd::Zero(tempData.size(), tempData[0].size());
        for (size_t i = 0; i < tempData.size(); ++i) {
            for (size_t j = 0; j < tempData[i].size(); ++j) {
                data.expressionMatrix(i, j) = tempData[i][j];
            }
        }

        // Load labels
        std::ifstream labFile(labelsFile);
        if (!labFile.is_open()) {
            throw DataLoadingError("Cannot open labels file: " + labelsFile);
        }

        std::getline(labFile, line); // Skip header
        int idx = 0;
        while (std::getline(labFile, line)) {
            std::istringstream iss(line);
            std::string sample, label;
            std::getline(iss, sample, ',');
            std::getline(iss, label, ',');
            data.labelIndices[label].push_back(idx++);
        }

        // Update metadata
        data.numGenes = data.expressionMatrix.rows();
        data.numSamples = data.expressionMatrix.cols();

        isDataLoaded = true;
        return true;
    }
    catch (const std::exception& e) {
        addError("Error loading data: " + std::string(e.what()));
        return false;
    }
}

void GeneExpressionAnalysis::preprocessData(const AnalysisParameters& params) {
    if (!isDataLoaded) {
        throw ProcessingError("Data not loaded");
    }

    validateParameters(params);

    try {
        // Log transform if requested and not already done
        if (params.performLogTransform && !state.isLogTransformed) {
            logTransform();
            state.isLogTransformed = true;
        }

        // Filter low expression genes
        if (!state.isFiltered) {
            filterLowExpression(params.minMeanExpression);
            filterLowVariance(params.minVariance);
            state.isFiltered = true;
        }

        // Normalize data
        if (params.performQuantileNorm && !state.isNormalized) {
            normalizeData();
            state.isNormalized = true;
        }

        // Center the data
        centerData();

        // Calculate QC metrics
        calculateQCMetrics();

        isPreprocessed = true;
        state.lastUsedParams = params;

    } catch (const std::exception& e) {
        addError("Error in preprocessing: " + std::string(e.what()));
        throw;
    }
}

void GeneExpressionAnalysis::logTransform() {
    // Add small constant to avoid log(0)
    double epsilon = 1e-10;
    data.expressionMatrix = (data.expressionMatrix.array() + epsilon).log2();
}

void GeneExpressionAnalysis::normalizeData() {
    // Quantile normalization
    int nRows = data.expressionMatrix.rows();
    int nCols = data.expressionMatrix.cols();

    // Sort each column
    Eigen::MatrixXd sortedMatrix = data.expressionMatrix;
    for (int j = 0; j < nCols; ++j) {
        std::vector<double> col(nRows);
        for (int i = 0; i < nRows; ++i) {
            col[i] = sortedMatrix(i, j);
        }
        std::sort(col.begin(), col.end());
        for (int i = 0; i < nRows; ++i) {
            sortedMatrix(i, j) = col[i];
        }
    }

    // Calculate row means
    Eigen::VectorXd rowMeans = sortedMatrix.rowwise().mean();

    // Create rank matrix
    Eigen::MatrixXd rankMatrix = data.expressionMatrix;
    for (int j = 0; j < nCols; ++j) {
        std::vector<std::pair<double, int>> colWithIdx(nRows);
        for (int i = 0; i < nRows; ++i) {
            colWithIdx[i] = {data.expressionMatrix(i, j), i};
        }
        std::sort(colWithIdx.begin(), colWithIdx.end());
        for (int i = 0; i < nRows; ++i) {
            rankMatrix(colWithIdx[i].second, j) = rowMeans(i);
        }
    }

    data.expressionMatrix = rankMatrix;
}

void GeneExpressionAnalysis::filterLowExpression(double threshold) {
    Eigen::VectorXd rowMeans = data.expressionMatrix.rowwise().mean();
    std::vector<int> validIndices;
    
    for (int i = 0; i < rowMeans.size(); ++i) {
        if (rowMeans(i) > threshold) {
            validIndices.push_back(i);
        }
    }

    // Create filtered matrix and update gene names
    Eigen::MatrixXd filteredMatrix(validIndices.size(), data.expressionMatrix.cols());
    std::vector<std::string> filteredGenes;
    
    for (size_t i = 0; i < validIndices.size(); ++i) {
        filteredMatrix.row(i) = data.expressionMatrix.row(validIndices[i]);
        filteredGenes.push_back(data.geneNames[validIndices[i]]);
    }

    data.expressionMatrix = filteredMatrix;
    data.geneNames = filteredGenes;
    data.numLowExpressionGenes = rowMeans.size() - validIndices.size();
}

void GeneExpressionAnalysis::centerData() {
    for (int i = 0; i < data.expressionMatrix.rows(); ++i) {
        double mean = data.expressionMatrix.row(i).mean();
        data.expressionMatrix.row(i).array() -= mean;
    }
}

Eigen::MatrixXd GeneExpressionAnalysis::calculateCorrelationMatrix(
    AnalysisParameters::CorrelationType type) const {
    
    if (!isPreprocessed) {
        throw ProcessingError("Data not preprocessed");
    }

    switch (type) {
        case AnalysisParameters::CorrelationType::PEARSON:
            return calculatePearsonCorrelation();
        case AnalysisParameters::CorrelationType::SPEARMAN:
            return calculateSpearmanCorrelation();
        default:
            throw ProcessingError("Unknown correlation type");
    }
}

Eigen::MatrixXd GeneExpressionAnalysis::calculatePearsonCorrelation() const {
    int numGenes = data.expressionMatrix.rows();
    Eigen::MatrixXd correlationMatrix = Eigen::MatrixXd::Zero(numGenes, numGenes);

    // Calculate correlations using matrix multiplication
    correlationMatrix = (data.expressionMatrix * data.expressionMatrix.transpose()) 
                       / (data.expressionMatrix.cols() - 1);

    return correlationMatrix;
}

void GeneExpressionAnalysis::calculateQCMetrics() {
    data.meanExpression = data.expressionMatrix.mean();
    data.medianExpression = [this]() {
        std::vector<double> values(data.expressionMatrix.data(), 
                                 data.expressionMatrix.data() + data.expressionMatrix.size());
        std::sort(values.begin(), values.end());
        return values[values.size() / 2];
    }();
}

void GeneExpressionAnalysis::exportResults(const ClusteringResult& results,
                                         const std::string& outputFile) const {
    std::ofstream out(outputFile);
    if (!out.is_open()) {
        throw ProcessingError("Cannot open output file: " + outputFile);
    }

    out << "Cluster,Gene,Mean Expression,Cluster Score\n";
    
    for (size_t i = 0; i < results.clusters.size(); ++i) {
        const auto& cluster = results.clusters[i];
        double clusterScore = results.clusterScores[i];
        
        for (int geneIdx : cluster) {
            if (geneIdx >= 0 && geneIdx < data.geneNames.size()) {
                double meanExp = data.expressionMatrix.row(geneIdx).mean();
                out << i << ","
                    << data.geneNames[geneIdx] << ","
                    << meanExp << ","
                    << clusterScore << "\n";
            }
        }
    }
}

void GeneExpressionAnalysis::exportCorrelationMatrix(const std::string& outputFile) const {
    std::ofstream out(outputFile);
    if (!out.is_open()) {
        throw ProcessingError("Cannot open output file: " + outputFile);
    }

    // Write header
    out << "Gene";
    for (const auto& gene : data.geneNames) {
        out << "," << gene;
    }
    out << "\n";

    // Write correlation values
    for (int i = 0; i < data.expressionMatrix.rows(); ++i) {
        out << data.geneNames[i];
        for (int j = 0; j < data.expressionMatrix.rows(); ++j) {
            out << "," << data.expressionMatrix(i, j);
        }
        out << "\n";
    }
}

void GeneExpressionAnalysis::validateParameters(const AnalysisParameters& params) const {
    if (params.minMeanExpression < 0) {
        throw std::invalid_argument("Minimum mean expression must be non-negative");
    }
    if (params.minVariance < 0) {
        throw std::invalid_argument("Minimum variance must be non-negative");
    }
    if (params.correlationThreshold < -1 || params.correlationThreshold > 1) {
        throw std::invalid_argument("Correlation threshold must be between -1 and 1");
    }
    if (params.minClusterSize < 1) {
        throw std::invalid_argument("Minimum cluster size must be positive");
    }
    if (params.maxClusterSize < params.minClusterSize) {
        throw std::invalid_argument("Maximum cluster size must be greater than minimum cluster size");
    }
}

} // namespace geneexp