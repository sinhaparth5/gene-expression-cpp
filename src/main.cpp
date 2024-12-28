#include "GeneExpressionAnalysis.hpp"
#include "ClusteringMethods.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <filesystem>

using namespace geneexp;

void printProgress(const std::string& stage, int progress, int total) {
    int barWidth = 50;
    float ratio = static_cast<float>(progress) / total;
    int filled = static_cast<int>(barWidth * ratio);
    
    std::cout << "\r" << stage << " [";
    for (int i = 0; i < barWidth; ++i) {
        if (i < filled) std::cout << "=";
        else std::cout << " ";
    }
    std::cout << "] " << int(ratio * 100.0) << "%" << std::flush;
}

int main(int argc, char* argv[]) {
    try {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " <expression_file> <labels_file>\n";
            return 1;
        }

        std::string expressionFile = argv[1];
        std::string labelsFile = argv[2];
        
        // Create output directory if it doesn't exist
        std::filesystem::create_directory("output");

        // Initialize analysis
        std::cout << "Initializing gene expression analysis...\n";
        GeneExpressionAnalysis analyzer;

        // Load data
        std::cout << "Loading data...\n";
        if (!analyzer.loadData(expressionFile, labelsFile)) {
            std::cerr << "Failed to load data\n";
            return 1;
        }

        // Set analysis parameters
        AnalysisParameters params;
        params.minMeanExpression = 2.0;
        params.minVariance = 0.2;
        params.performLogTransform = true;
        params.performQuantileNorm = true;
        params.correlationType = AnalysisParameters::CorrelationType::PEARSON;
        params.correlationThreshold = 0.4;
        params.minClusterSize = 10;
        params.maxClusterSize = 500;

        // Preprocess data
        std::cout << "Preprocessing data...\n";
        analyzer.preprocessData(params);

        // Calculate correlation matrix
        std::cout << "Calculating correlation matrix...\n";
        Eigen::MatrixXd correlationMatrix = analyzer.calculateCorrelationMatrix(params.correlationType);

        // Perform clustering
        std::cout << "Performing clustering analysis...\n";
        auto clusteringMethod = clustering::ClusteringMethodFactory::createMethod("hierarchical", params);
        ClusteringResult results = clusteringMethod->performClustering(correlationMatrix, params);

        // Export results
        std::cout << "Exporting results...\n";
        analyzer.exportResults(results, "output/gene_clusters.csv");
        analyzer.exportCorrelationMatrix("output/correlation_matrix.csv");
        analyzer.exportClusterVisualization(results, "output/cluster_visualization.html");

        // Print summary
        std::cout << "\nAnalysis complete!\n";
        std::cout << "Number of clusters: " << results.clusters.size() << "\n";
        std::cout << "Results saved in the 'output' directory\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}