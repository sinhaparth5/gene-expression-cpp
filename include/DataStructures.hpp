#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>
#include <memory>

namespace geneexp {

// Structure to hold row gene expression data
struct GeneExpressionData {
    Eigen::MatrixXd expressionMatrix; // Rows are genes, columns are samples
    std::vector<std::string> geneNames;
    std::vector<std::string> sampleLabels;
    std::map<std::string, std::vector<int>> labelIndices;   // Maps class labels to sample indices

    // Metadata
    std::string datasetName;
    std::string description;
    int numGenes;
    int numSamples;

    // Quality control metrics
    double medianExpression;
    double meanExpression;
    int numLowExpressionGenes;
    int numHighVarianceGenes;
};

// Parameters for data preprocessing and analysis
struct AnalysisParameters {
    // Preproccessing parameters
    double minMeanExpression = 1.0;
    double minVariance = 0.1;
    bool performLogTransform = true;
    bool performQuantileNorm = true;

    //Correlation parameters
    enum class CorrelationType {
        PEARSON,
        SPEARMAN
    } correlationType = CorrelationType::PEARSON;

    // Clustering parameters
    double correlationThreshold = 0.3;
    int minClusterSize = 5;
    int maxClusterSize = 1000;
    double mergingThreshold = 0.8;
};

//Structure to hold clustering results
struct ClusteringResult {
    std::vector<std::vector<int>> clusters;
    std::vector<double> clusterScores;
    Eigen::MatrixXd correlationMatrix;

    // Cluster statistics
    struct ClusterStats {
        double meanCorrelation;
        double minCorrelation;
        double maxCorrelation;
        std::vector<std::string> enrichedPathways;
        std::vector<double> pathwayPValues;
    };
    std::vector<ClusterStats> clusterStatistics;
};

// Structure for starting analysis results
struct AnalysisResult {
    GeneExpressionData processedData;
    ClusteringResult clusteringResult;
    AnalysisParameters usedParameters;

    // Analysis metadata
    std::string timestamp;
    std::string analysisVersion;
    std::vector<std::string> warnings;
    std::vector<std::string> errors;
};

}