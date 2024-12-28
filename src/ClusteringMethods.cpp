#include "ClusteringMethods.hpp"
#include <algorithm>
#include <limits>
#include <random>
#include <cmath>
#include <numeric>

namespace geneexp {
namespace clustering {

// ... [previous code remains the same]

double HierarchicalClustering::calculateLinkageDistance(
    const std::vector<double>& distances,
    LinkageType type) {
    
    switch (type) {
        case LinkageType::SINGLE:
            return *std::min_element(distances.begin(), distances.end());
            
        case LinkageType::COMPLETE:
            return *std::max_element(distances.begin(), distances.end());
            
        case LinkageType::AVERAGE: {
            double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
            return sum / distances.size();
        }
            
        case LinkageType::WARD: {
            double sum = 0.0;
            for (double d : distances) {
                sum += d * d;
            }
            return std::sqrt(sum / distances.size());
        }
            
        default:
            throw std::runtime_error("Unknown linkage type");
    }
}

std::string HierarchicalClustering::getDescription() const {
    std::string desc = "Hierarchical clustering with ";
    switch (linkageType) {
        case LinkageType::SINGLE:
            desc += "single linkage";
            break;
        case LinkageType::COMPLETE:
            desc += "complete linkage";
            break;
        case LinkageType::AVERAGE:
            desc += "average linkage";
            break;
        case LinkageType::WARD:
            desc += "Ward's method";
            break;
    }
    return desc;
}

// K-means Clustering Implementation
KMeansClustering::KMeansClustering(int k, int maxIterations)
    : k(k), maxIterations(maxIterations) {}

ClusteringResult KMeansClustering::performClustering(
    const Eigen::MatrixXd& correlationMatrix,
    const AnalysisParameters& params) {
    
    int n = correlationMatrix.rows();
    
    // Initialize centroids
    Eigen::MatrixXd centroids = initializeCentroids(correlationMatrix, k);
    std::vector<int> assignments(n);
    
    // Iterate until convergence or max iterations
    for (int iter = 0; iter < maxIterations; ++iter) {
        // Assign points to nearest centroid
        std::vector<int> newAssignments = assignClusters(correlationMatrix, centroids);
        
        // Check for convergence
        if (newAssignments == assignments) {
            break;
        }
        
        assignments = newAssignments;
        
        // Update centroids
        centroids = updateCentroids(correlationMatrix, assignments, k);
    }
    
    // Convert assignments to cluster format
    std::vector<std::vector<int>> clusters(k);
    for (int i = 0; i < n; ++i) {
        clusters[assignments[i]].push_back(i);
    }
    
    // Calculate cluster scores
    std::vector<double> clusterScores;
    std::vector<ClusteringResult::ClusterStats> clusterStats;
    
    for (const auto& cluster : clusters) {
        if (cluster.size() >= params.minClusterSize) {
            double totalCorr = 0.0;
            double minCorr = 1.0;
            double maxCorr = -1.0;
            int pairs = 0;
            
            for (size_t i = 0; i < cluster.size(); ++i) {
                for (size_t j = i + 1; j < cluster.size(); ++j) {
                    double corr = correlationMatrix(cluster[i], cluster[j]);
                    totalCorr += corr;
                    minCorr = std::min(minCorr, corr);
                    maxCorr = std::max(maxCorr, corr);
                    pairs++;
                }
            }
            
            double meanCorr = pairs > 0 ? totalCorr / pairs : 0.0;
            clusterScores.push_back(meanCorr);
            
            ClusteringResult::ClusterStats stats{
                meanCorr,
                minCorr,
                maxCorr,
                std::vector<std::string>(),
                std::vector<double>()
            };
            clusterStats.push_back(stats);
        }
    }
    
    return ClusteringResult{
        clusters,
        clusterScores,
        correlationMatrix,
        clusterStats
    };
}

Eigen::MatrixXd KMeansClustering::initializeCentroids(
    const Eigen::MatrixXd& data,
    int k) {
    
    // Use k-means++ initialization
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    int n = data.rows();
    std::vector<int> centroidIndices;
    centroidIndices.push_back(std::floor(dis(gen) * n));
    
    std::vector<double> distances(n, std::numeric_limits<double>::infinity());
    
    // Choose remaining centroids
    for (int i = 1; i < k; ++i) {
        // Update distances
        for (int j = 0; j < n; ++j) {
            double minDist = std::numeric_limits<double>::infinity();
            for (int c : centroidIndices) {
                double dist = (data.row(j) - data.row(c)).squaredNorm();
                minDist = std::min(minDist, dist);
            }
            distances[j] = minDist;
        }
        
        // Choose next centroid
        double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
        double r = dis(gen) * sum;
        double cumSum = 0.0;
        int nextCentroid = 0;
        
        for (int j = 0; j < n; ++j) {
            cumSum += distances[j];
            if (cumSum >= r) {
                nextCentroid = j;
                break;
            }
        }
        
        centroidIndices.push_back(nextCentroid);
    }
    
    // Create centroid matrix
    Eigen::MatrixXd centroids(k, data.cols());
    for (int i = 0; i < k; ++i) {
        centroids.row(i) = data.row(centroidIndices[i]);
    }
    
    return centroids;
}

std::vector<int> KMeansClustering::assignClusters(
    const Eigen::MatrixXd& data,
    const Eigen::MatrixXd& centroids) {
    
    int n = data.rows();
    std::vector<int> assignments(n);
    
    for (int i = 0; i < n; ++i) {
        double minDist = std::numeric_limits<double>::infinity();
        int bestCluster = 0;
        
        for (int j = 0; j < centroids.rows(); ++j) {
            double dist = (data.row(i) - centroids.row(j)).squaredNorm();
            if (dist < minDist) {
                minDist = dist;
                bestCluster = j;
            }
        }
        
        assignments[i] = bestCluster;
    }
    
    return assignments;
}

Eigen::MatrixXd KMeansClustering::updateCentroids(
    const Eigen::MatrixXd& data,
    const std::vector<int>& assignments,
    int k) {
    
    Eigen::MatrixXd newCentroids = Eigen::MatrixXd::Zero(k, data.cols());
    std::vector<int> counts(k, 0);
    
    // Sum points in each cluster
    for (size_t i = 0; i < assignments.size(); ++i) {
        int cluster = assignments[i];
        newCentroids.row(cluster) += data.row(i);
        counts[cluster]++;
    }
    
    // Calculate means
    for (int i = 0; i < k; ++i) {
        if (counts[i] > 0) {
            newCentroids.row(i) /= counts[i];
        }
    }
    
    return newCentroids;
}

std::string KMeansClustering::getDescription() const {
    return "K-means clustering with k=" + std::to_string(k) + 
           " and max_iterations=" + std::to_string(maxIterations);
}

// Factory Implementation
std::unique_ptr<ClusteringMethod> ClusteringMethodFactory::createMethod(
    const std::string& methodName,
    const AnalysisParameters& params) {
    
    if (methodName == "hierarchical") {
        return std::make_unique<HierarchicalClustering>(
            HierarchicalClustering::LinkageType::COMPLETE);
    } else if (methodName == "kmeans") {
        return std::make_unique<KMeansClustering>();
    } else {
        throw std::runtime_error("Unknown clustering method: " + methodName);
    }
}

std::vector<std::string> ClusteringMethodFactory::availableMethods() {
    return {"hierarchical", "kmeans"};
}

} // namespace clustering
} // namespace geneexp