#include "meshedparticle.h"
#include "meshops.h"
#include "energy.h"

void ParticleAdhesion::find_pairs(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd V_particle, Eigen::MatrixXi F_particle,double distance_threshold, std::vector<std::pair<int, int>> &bonds)
{
Mesh M1;
 //const double distance_threshold = 0.20;
nearest.resize(V.rows());
nearest.setConstant(-1); // initialize to -1 to mark vertices that don't have a nearest neighbor
for (int i = 0; i < V.rows(); i++) {
    double min_distance = std::numeric_limits<double>::max();
    int nearest_index = -1;
    for (int j = 0; j < V_particle.rows(); j++) {
        double distance = sqrt(pow(V(i,0)-V_particle(j,0), 2)
                             +pow(V(i,1)-V_particle(j,1), 2)
                             +pow(V(i,2)-V_particle(j,2), 2));
        if (distance < min_distance) {
            min_distance = distance;
            nearest_index = j;
        }
    }
    if (nearest_index != -1 && min_distance < distance_threshold) {
        nearest(i) = nearest_index;
    }
}
//std::cout << nearest<< std::endl;
// Create bond connections based on nearest neighbors
//std::vector<std::pair<int, int>> bonds;
for (int i = 0; i < V.rows(); i++) {
    int j = nearest(i);
    if (j != -1) {
        bonds.emplace_back(i, j);
    }
}

// Output the number of bonds and the list of bonds
// std::cout << "Number of bonds: " << bonds.size() << std::endl;
// for (const auto& bond : bonds) {
//     std::cout << bond.first << " " << bond.second << std::endl;
// }   
//return bonds;
}

/*
#include <igl/kdtree/kdtree.h>

void find_pairs_2(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2, double distance_threshold, std::vector<std::pair<int, int>>& bonds)
{
    // Build a k-d tree for V2
    igl::AABB<Eigen::MatrixXd, 3> tree;
    tree.init(V2);

    // Search for nearest neighbors of each vertex in V1 in the k-d tree
    std::vector<int> nearest_indices;
    std::vector<double> nearest_dists2;
    nearest_indices.resize(V1.rows());
    nearest_dists2.resize(V1.rows());
    tree.squared_distance(V1, nearest_indices, nearest_dists2);

    // Create bond connections based on nearest neighbors
    bonds.clear();
    for (int i = 0; i < V1.rows(); i++) {
        if (nearest_dists2[i] < distance_threshold*distance_threshold) {
            bonds.emplace_back(i, nearest_indices[i]);
        }
    }
}
*/