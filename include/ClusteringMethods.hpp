#pragma once

#include "DataStructures.hpp"
#include <Eigen/Dense>
#include <vector>
#include <memory>

namespace geneexp {
namespace clustering {

class ClusteringMethod {
public:
    virtual ~ClusteringMethod() = default;
    
    virtual ClusteringResult performClustering(
        const Eigen::MatrixXd& correlationMatrix,
        const AnalysisParameters& params) = 0;
        
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
};

class HierarchicalClustering : public ClusteringMethod {
public:
    enum class LinkageType {
        SINGLE,
        COMPLETE,
        AVERAGE,
        WARD
    };

    explicit HierarchicalClustering(LinkageType linkage = LinkageType::COMPLETE);
    
    ClusteringResult performClustering(
        const Eigen::MatrixXd& correlationMatrix,
        const AnalysisParameters& params) override;
        
    std::string getName() const override { return "Hierarchical Clustering"; }
    std::string getDescription() const override;

private:
    LinkageType linkageType;
    
    std::vector<std::vector<int>> mergeClusters(
        const std::vector<std::vector<int>>& clusters,
        const Eigen::MatrixXd& distanceMatrix,
        double threshold);
        
    double calculateClusterDistance(
        const std::vector<int>& cluster1,
        const std::vector<int>& cluster2,
        const Eigen::MatrixXd& distanceMatrix);
        
    double calculateLinkageDistance(
        const std::vector<double>& distances,
        LinkageType type);
};

class KMeansClustering : public ClusteringMethod {
public:
    explicit KMeansClustering(int k = 10, int maxIter = 100);
    
    ClusteringResult performClustering(
        const Eigen::MatrixXd& correlationMatrix,
        const AnalysisParameters& params) override;
        
    std::string getName() const override { return "K-means Clustering"; }
    std::string getDescription() const override;

private:
    int k;
    int maxIter;  // Changed from maxIterations
    
    Eigen::MatrixXd initializeCentroids(
        const Eigen::MatrixXd& data,
        int k);
        
    std::vector<int> assignClusters(
        const Eigen::MatrixXd& data,
        const Eigen::MatrixXd& centroids);
        
    Eigen::MatrixXd updateCentroids(
        const Eigen::MatrixXd& data,
        const std::vector<int>& assignments,
        int k);
};

class ClusteringMethodFactory {
public:
    static std::unique_ptr<ClusteringMethod> createMethod(
        const std::string& methodName,
        const AnalysisParameters& params);
        
    static std::vector<std::string> availableMethods();
};

} // namespace clustering
} // namespace geneexp