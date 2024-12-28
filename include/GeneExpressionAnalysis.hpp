#pragma once

#include "DataStructures.hpp"
#include <string>
#include <memory>
#include <stdexcept>

namespace geneexp {

class GeneExpressionAnalysis {
public:
    GeneExpressionAnalysis();
    ~GeneExpressionAnalysis();

    // Data loading and validation
    bool loadData(const std::string& expressionFile, const std::string& labelsFile);
    bool validateData() const;

    // Preprocessing pipeline
    void preprocessData(const AnalysisParameters& param);

    // Individual preprocessing steps
    void normalizeData();
    void filterLowExpression(double threshold);
    void filterLowVariance(double threshold);
    void logTransform();
    void quantileNormalization();

    //Analysis methods
    Eigen::MatrixXd calculateCorrelationMatrix(AnalysisParameters::CorrelationType type) const;
    ClusteringResult performClustering(const AnalysisParameters& params);

    // Statistical analysis
    void calculateClusterStatistics(ClusteringResult& result) const;
    void performPathwayEnrichment(ClusteringResult& result) const;

    // Quality
    void calculateQCMetrics();
    bool checkQualityControl() const;

    // Export methods
    void exportResults(const ClusteringResult& results, const std::string& outputDir) const;
    void exportCorrelationMatrix(const std::string& outputFile) const;
    void exportClusterVisualization(const ClusteringResult& results, const std::string& outputFile) const;

    // Add the new visualization method with the parameters
    void exportClusterVisualization(
        const ClusteringResult& results,
        const std::string& outputFile,
        int maxGenesPerCluster = 3
    ) const;

    // Image export
    void exportBestGeneImage(const std::string& outputFile) const;

    // Getters
    const GeneExpressionData& getData() const { return data; }
    const std::vector<std::string>& getWarnings() const { return warnings; }
    const std::vector<std::string>& getErrors() const { return errors; }

private:
    GeneExpressionData data;
    std::vector<std::string> warnings;
    std::vector<std::string> errors;
    bool isDataLoaded;
    bool isPreprocessed;
    std::vector<int> selectTopGenes(int numGenes) const;

    // Helper methods
    void centerData();
    double calculateGeneVariance(const Eigen::VectorXd& geneExpression) const;
    Eigen::MatrixXd calculateSpearmanCorrelation() const;
    Eigen::MatrixXd calculatePearsonCorrelation() const;
    void validateParameters(const AnalysisParameters& params) const;
    void addWarning(const std::string& warning) { warnings.push_back(warning); }
    void addError(const std::string& error) { errors.push_back(error); }

    // Error handling
    struct ProcessingState {
        bool isNormalized;
        bool isLogTransformed;
        bool isFiltered;
        AnalysisParameters lastUsedParams;
    } state;
};

// Custom exceptions
class GeneExpressionError : public std::runtime_error {
public:
    explicit GeneExpressionError(const std::string& message)
        : std::runtime_error(message) {}
};

class DataLoadingError : public GeneExpressionError {
public:
    explicit DataLoadingError(const std::string& message)
        : GeneExpressionError(message) {}
};

class ProcessingError : public GeneExpressionError {
public:
    explicit ProcessingError(const std::string& message)
        : GeneExpressionError(message) {}
};


}