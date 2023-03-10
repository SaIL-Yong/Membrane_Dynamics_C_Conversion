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
std::cout << nearest<< std::endl;
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