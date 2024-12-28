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
        std::cout << "Starting to load data files..." << std::endl;
        
        // Load expression data
        std::ifstream expFile(expressionFile);
        if (!expFile.is_open()) {
            std::cerr << "Error: Cannot open expression file: " << expressionFile << std::endl;
            return false;
        }

        // Read header
        std::string line;
        if (!std::getline(expFile, line)) {
            std::cerr << "Error: Expression file is empty" << std::endl;
            return false;
        }

        std::istringstream iss(line);
        std::string token;
        std::getline(iss, token, ','); // Skip the empty first column header
        
        // Parse header to get gene names (columns are genes in this case)
        while (std::getline(iss, token, ',')) {
            data.geneNames.push_back(token);
        }

        if (data.geneNames.empty()) {
            std::cerr << "Error: No genes found in header" << std::endl;
            return false;
        }
        std::cout << "Found " << data.geneNames.size() << " genes in header" << std::endl;

        // Read sample data
        std::vector<std::vector<double>> tempData;
        int lineNum = 1;
        while (std::getline(expFile, line)) {
            lineNum++;
            std::vector<double> row;
            std::istringstream iss(line);
            
            // Get sample name
            std::string sampleName;
            if (!std::getline(iss, sampleName, ',')) {
                std::cerr << "Error: Invalid format at line " << lineNum << std::endl;
                continue;
            }
            data.sampleLabels.push_back(sampleName);
            
            // Get expression values
            while (std::getline(iss, token, ',')) {
                try {
                    double value = std::stod(token);
                    row.push_back(value);
                } catch (const std::exception& e) {
                    std::cerr << "Warning: Invalid value '" << token << "' at line " << lineNum 
                             << ", using 0.0" << std::endl;
                    row.push_back(0.0);
                }
            }

            if (row.size() != data.geneNames.size()) {
                std::cerr << "Error: Inconsistent number of values at line " << lineNum 
                         << ". Expected: " << data.geneNames.size() 
                         << ", Got: " << row.size() << std::endl;
                continue;
            }

            tempData.push_back(row);
        }

        if (tempData.empty()) {
            std::cerr << "Error: No valid expression data found" << std::endl;
            return false;
        }

        // Convert to Eigen matrix (samples as rows, genes as columns)
        data.expressionMatrix = Eigen::MatrixXd::Zero(tempData.size(), data.geneNames.size());
        for (size_t i = 0; i < tempData.size(); ++i) {
            for (size_t j = 0; j < tempData[i].size(); ++j) {
                data.expressionMatrix(i, j) = tempData[i][j];
            }
        }

        std::cout << "Expression matrix dimensions: " 
                  << data.expressionMatrix.rows() << " samples × " 
                  << data.expressionMatrix.cols() << " genes" << std::endl;

        // Load labels
        std::ifstream labFile(labelsFile);
        if (!labFile.is_open()) {
            std::cerr << "Error: Cannot open labels file: " << labelsFile << std::endl;
            return false;
        }

        std::getline(labFile, line); // Skip header
        std::map<std::string, std::string> sampleToClass;
        
        while (std::getline(labFile, line)) {
            std::istringstream iss(line);
            std::string sample, label;
            
            std::getline(iss, sample, ','); // Skip empty first column
            if (!std::getline(iss, sample, ',')) continue;
            if (!std::getline(iss, label, ',')) continue;
            
            sampleToClass[sample] = label;
        }

        // Map samples to their class labels
        for (size_t i = 0; i < data.sampleLabels.size(); ++i) {
            std::string label = sampleToClass[data.sampleLabels[i]];
            data.labelIndices[label].push_back(i);
        }

        // Update metadata
        data.numGenes = data.expressionMatrix.cols();
        data.numSamples = data.expressionMatrix.rows();
        data.datasetName = "Gene Expression Dataset";
        data.description = "Dataset with " + std::to_string(data.numGenes) + 
                          " genes across " + std::to_string(data.numSamples) + " samples";

        std::cout << "Successfully loaded data with:" << std::endl
                  << "- " << data.numSamples << " samples" << std::endl
                  << "- " << data.numGenes << " genes" << std::endl
                  << "- " << data.labelIndices.size() << " unique classes" << std::endl;

        isDataLoaded = true;
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "Error loading data: " << e.what() << std::endl;
        return false;
    }
}

void GeneExpressionAnalysis::preprocessData(const AnalysisParameters& params) {
    if (!isDataLoaded) {
        throw ProcessingError("Data not loaded");
    }

    validateParameters(params);

    try {
        // Add progress reporting
        std::cout << "Starting preprocessing..." << std::endl;
        
        // Pre-filter genes with very low expression to reduce memory usage early
        if (!state.isFiltered) {
            std::cout << "Filtering low expression genes..." << std::endl;
            filterLowExpression(params.minMeanExpression);
            
            // Only proceed with variance filtering if we have enough genes left
            if (data.numGenes > 100) {
                std::cout << "Filtering low variance genes..." << std::endl;
                filterLowVariance(params.minVariance);
            }
            state.isFiltered = true;
        }

        // Log transform if requested and not already done
        if (params.performLogTransform && !state.isLogTransformed) {
            std::cout << "Performing log transformation..." << std::endl;
            logTransform();
            state.isLogTransformed = true;
        }

        // Normalize data
        if (params.performQuantileNorm && !state.isNormalized) {
            std::cout << "Performing normalization..." << std::endl;
            normalizeData();
            state.isNormalized = true;
        }

        // Center the data
        std::cout << "Centering data..." << std::endl;
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
    // Calculate column-wise means since genes are columns
    Eigen::VectorXd geneMeans = data.expressionMatrix.colwise().mean();
    std::vector<int> validIndices;
    
    for (int i = 0; i < geneMeans.size(); ++i) {
        if (geneMeans(i) > threshold) {
            validIndices.push_back(i);
        }
    }

    // Create filtered matrix maintaining sample rows × gene columns structure
    Eigen::MatrixXd filteredMatrix(data.expressionMatrix.rows(), validIndices.size());
    std::vector<std::string> filteredGenes;
    
    for (size_t i = 0; i < validIndices.size(); ++i) {
        filteredMatrix.col(i) = data.expressionMatrix.col(validIndices[i]);
        filteredGenes.push_back(data.geneNames[validIndices[i]]);
    }

    data.expressionMatrix = filteredMatrix;
    data.geneNames = filteredGenes;
    data.numGenes = validIndices.size();
    data.numLowExpressionGenes = geneMeans.size() - validIndices.size();
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

void GeneExpressionAnalysis::filterLowVariance(double threshold) {
    std::vector<int> validIndices;
    validIndices.reserve(data.numGenes);

    // Calculate variance for each gene (column)
    for (int i = 0; i < data.expressionMatrix.cols(); ++i) {
        Eigen::VectorXd col = data.expressionMatrix.col(i);
        double mean = col.mean();
        double variance = (col.array() - mean).square().sum() / (col.size() - 1);
        
        if (variance > threshold) {
            validIndices.push_back(i);
        }
    }

    // Create filtered matrix maintaining sample rows × gene columns structure
    Eigen::MatrixXd filteredMatrix(data.expressionMatrix.rows(), validIndices.size());
    std::vector<std::string> filteredGenes;
    filteredGenes.reserve(validIndices.size());

    for (size_t i = 0; i < validIndices.size(); ++i) {
        filteredMatrix.col(i) = data.expressionMatrix.col(validIndices[i]);
        filteredGenes.push_back(data.geneNames[validIndices[i]]);
    }

    data.numHighVarianceGenes = validIndices.size();
    data.expressionMatrix = filteredMatrix;
    data.geneNames = filteredGenes;
    data.numGenes = validIndices.size();
}

Eigen::MatrixXd GeneExpressionAnalysis::calculateSpearmanCorrelation() const {
    const int numGenes = data.expressionMatrix.rows();
    const int numSamples = data.expressionMatrix.cols();
    Eigen::MatrixXd rankMatrix(numGenes, numSamples);

    // Calculate ranks for each gene
    for (int i = 0; i < numGenes; ++i) {
        std::vector<std::pair<double, int>> values(numSamples);
        for (int j = 0; j < numSamples; ++j) {
            values[j] = {data.expressionMatrix(i, j), j};
        }
        
        std::sort(values.begin(), values.end());
        
        // Handle ties by using average ranks
        for (int j = 0; j < numSamples;) {
            int k = j + 1;
            while (k < numSamples && values[k].first == values[j].first) {
                k++;
            }
            
            // Calculate average rank for ties
            double avgRank = (j + k - 1) / 2.0 + 1;
            for (int l = j; l < k; l++) {
                rankMatrix(i, values[l].second) = avgRank;
            }
            j = k;
        }
    }

    // Calculate correlation on ranks using Pearson correlation formula
    Eigen::MatrixXd correlationMatrix = Eigen::MatrixXd::Zero(numGenes, numGenes);
    
    for (int i = 0; i < numGenes; ++i) {
        Eigen::VectorXd ranksI = rankMatrix.row(i);
        double meanI = ranksI.mean();
        double stdI = std::sqrt((ranksI.array() - meanI).square().sum() / (numSamples - 1));

        for (int j = i; j < numGenes; ++j) {
            Eigen::VectorXd ranksJ = rankMatrix.row(j);
            double meanJ = ranksJ.mean();
            double stdJ = std::sqrt((ranksJ.array() - meanJ).square().sum() / (numSamples - 1));

            if (stdI > 0 && stdJ > 0) {
                double correlation = ((ranksI.array() - meanI) * (ranksJ.array() - meanJ)).sum() 
                                   / ((numSamples - 1) * stdI * stdJ);
                correlationMatrix(i, j) = correlation;
                correlationMatrix(j, i) = correlation;
            }
        }
    }

    return correlationMatrix;
}

std::vector<int> GeneExpressionAnalysis::selectTopGenes(int numGenes) const {
    std::vector<std::pair<double, int>> geneVariances;
    
    // Calculate variance for each gene
    for (int i = 0; i < data.expressionMatrix.cols(); ++i) {
        Eigen::VectorXd col = data.expressionMatrix.col(i);
        double mean = col.mean();
        double variance = (col.array() - mean).square().sum() / (col.size() - 1);
        geneVariances.push_back({variance, i});
    }
    
    // Sort by variance in descending order
    std::sort(geneVariances.begin(), geneVariances.end(),
              std::greater<std::pair<double, int>>());
    
    // Select top N genes
    std::vector<int> topGenes;
    for (int i = 0; i < std::min(numGenes, (int)geneVariances.size()); ++i) {
        topGenes.push_back(geneVariances[i].second);
    }
    
    return topGenes;
}

void GeneExpressionAnalysis::exportBestGeneImage(const std::string& outputFile) const {
    const int IMAGE_WIDTH = 1000;
    const int IMAGE_HEIGHT = 200;

    std::vector<unsigned char> image(IMAGE_WIDTH * IMAGE_HEIGHT * 3, 255);

    // Find the gene with highest variance (reuse same logic)
    double maxVariance = -1;
    int bestGeneIdx = 0;

    for (int i = 0; i < data.expressionMatrix.cols(); ++i) {
        Eigen::VectorXd col = data.expressionMatrix.col(i);
        double mean = col.mean();
        double variance = (col.array() - mean).square().sum() / (col.size() - 1);

        if (variance > maxVariance) {
            maxVariance = variance;
            bestGeneIdx = i;
        }
    }

    // Get expression values for best gene
    Eigen::VectorXd geneExpr = data.expressionMatrix.col(bestGeneIdx);

    // Calculate min and max for scaling
    double minExpr = geneExpr.minCoeff();
    double maxExpr = geneExpr.maxCoeff();

    double width = static_cast<double>(IMAGE_WIDTH) / geneExpr.size();

    for (int i = 0; i < geneExpr.size(); ++i) {
        //Normalize expression value to 0-1 range
        double normalizedExpr = (geneExpr(i) - minExpr) / (maxExpr - minExpr);

        // Convert to blue intensity
        int blueIntensity = static_cast<int>(255 * (1.0 - normalizedExpr));

        // Fill the column in the image
        int x = static_cast<int>(i * width);
        int nextX = static_cast<int>((i + 1) * width);

        for (int px = x; px < nextX && px < IMAGE_WIDTH; ++px) {
            for (int y = 50; y < 150; ++y) {
                int idx = (y * IMAGE_WIDTH + px) * 3;
                image[idx] = blueIntensity;
                image[idx + 1] = blueIntensity;
                image[idx + 2] = 255;
            }
        }
    }

    //Save as ppm (a simple image format)
    std::ofstream out(outputFile, std::ios::binary);
    if(!out.is_open()) {
        throw ProcessingError("Connot open output file: " + outputFile);
    }

    out << "P6\n" << IMAGE_WIDTH << " " << IMAGE_HEIGHT << "\n255\n";

    //Write image data
    out.write(reinterpret_cast<char*>(image.data()), image.size());
}

void GeneExpressionAnalysis::exportClusterVisualization(
    const ClusteringResult& results, 
    const std::string& outputFile,
    int maxGenesPerCluster) const {
    
    std::ofstream out(outputFile);
    if (!out.is_open()) {
        throw ProcessingError("Cannot open output file: " + outputFile);
    }

    // Write HTML header with simplified styling
    out << R"(<!DOCTYPE html>
<html>
<head>
    <title>Gene Expression Clusters (Subset)</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .heatmap { border-collapse: collapse; }
        .heatmap td { width: 3px; height: 25px; padding: 0; }
        .gene-label { 
            padding-right: 10px; 
            font-family: monospace; 
            background: white;
        }
        .cluster-header { 
            padding: 10px 0; 
            font-weight: bold; 
            background: #f0f0f0;
        }
    </style>
</head>
<body>)";

    out << "<h2>Gene Expression Clusters (Top Genes)</h2>"
        << "<p>Showing up to " << maxGenesPerCluster << " genes per cluster</p>";

    out << "<table class='heatmap'>\n";
    
    for (size_t i = 0; i < results.clusters.size(); ++i) {
        const auto& cluster = results.clusters[i];
        if (cluster.empty()) continue;

        // Get indices of genes in this cluster
        std::vector<std::pair<double, int>> clusterGeneVariances;
        for (int geneIdx : cluster) {
            Eigen::VectorXd col = data.expressionMatrix.col(geneIdx);
            double mean = col.mean();
            double variance = (col.array() - mean).square().sum() / (col.size() - 1);
            clusterGeneVariances.push_back({variance, geneIdx});
        }
        
        // Sort cluster genes by variance
        std::sort(clusterGeneVariances.begin(), clusterGeneVariances.end(),
                 std::greater<std::pair<double, int>>());
        
        // Take only top N genes from this cluster
        std::vector<int> clusterSubset;
        for (size_t j = 0; j < std::min(static_cast<size_t>(maxGenesPerCluster), 
                                      clusterGeneVariances.size()); ++j) {
            clusterSubset.push_back(clusterGeneVariances[j].second);
        }

        if (clusterSubset.empty()) continue;

        // Cluster header
        out << "<tr><td colspan='" << (data.expressionMatrix.cols() + 1) 
            << "' class='cluster-header'>Cluster " << i 
            << " (Showing " << clusterSubset.size() << " top genes)</td></tr>\n";
        
        // Output only the selected genes
        for (int geneIdx : clusterSubset) {
            out << "<tr>\n<td class='gene-label'>" << data.geneNames[geneIdx] << "</td>\n";
            
            // Add expression values
            for (int j = 0; j < data.expressionMatrix.cols(); ++j) {
                double value = data.expressionMatrix(geneIdx, j);
                int intensity = static_cast<int>((value + 4) * 32);
                intensity = std::min(255, std::max(0, intensity));
                
                out << "<td style='background-color: rgb(" 
                    << (255 - intensity) << "," 
                    << (255 - intensity) << "," 
                    << 255 << ")'></td>\n";
            }
            out << "</tr>\n";
        }
    }
    
    out << "</table>\n</body>\n</html>";
}
} // namespace geneexp