#pragma once

#include "DataStructures.hpp"
#include <Eigen/Dense>
#include <vector>
#include <memory>

namespace geneexp {
namespace clustering {

// Abstract base class for clustering methods
class ClusteringMethod {
public:
    virtual ~ClusteringMethod() = default;

    virtual ClusteringResult performClustering(
        const Eigen::MatrixXd& correlationMatrix,
        const AnalysisParameters& params
    ) = 0;

    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
};

// Hierarchical clustering implementation
class HierarchicalClustering : public ClusteringMethod {
public:
    enum class LinkageType {
        SINGLE,
        COMPLETE,
        AVERAGE,
        WARD
    };

    HierarchicalClustering(LinkageType linkage = LinkageType::COMPLETE);

    ClusteringResult performClustering(
        const Eigen::MatrixXd& correlationMatrix,
        const AnalysisParameters& params
    ) override;

    std::string getName() const override { return "Hierarchical Clustering"; }
    std::string getDescription() const override;

private:
    LinkageType linkageType;

    struct ClusterNode {
        std::vector<int> members;
        double height;
        std::shared_ptr<ClusterNode> left;
        std::shared_ptr<ClusterNode> right;
    };

    std::vector<std::vector<int>> mergeCluster(
        const std::vector<std::vector<int>>& clusters,
        const Eigen::MatrixXd& distanceMatrix,
        double threshold
    );

    double calculateClusterDistance(
        const std::vector<int>& cluster1,
        const std::vector<int>& cluster2,
        const Eigen::MatrixXd& distanceMatrix
    );

    double calculateLinkageDistance(
        const std::vector<double>& distance,
        LinkageType type
    );
};

//K-means clustering implementation
class KMeansClustering : public ClusteringMethod {
public:
    explicit KMeansClustering(int k = 10, int maxIterations = 100);

    ClusteringResult performClustering(
        const Eigen::MatrixXd& correlationMatrix,
        const AnalysisParameters& params
    ) override;

    std::string getName() const override { return "K-Means Clustering"; }
    std::string getDescription() const override;

private:
    int k;
    int maxIneration;

    Eigen::MatrixXd initializeCentroids(
        const Eigen::MatrixXd& data,
        int k
    );

    std::vector<int> assignCluster(
        const Eigen::MatrixXd& data,
        const Eigen::MatrixXd& centroids
    );

    Eigen::MatrixXd updateCentroids(
        const Eigen::MatrixXd& data,
        const std::vector<int>& assignments,
        int k
    );
};

class ClusteringMethodFactory {
public:
    static std::unique_ptr<ClusteringMethod> createMethod(
        const std::string& methodName,
        const AnalysisParameters& params
    );

    static std::vector<std::string> availableMethods();
};

}
}