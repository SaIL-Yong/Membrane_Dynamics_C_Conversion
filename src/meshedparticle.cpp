#include "meshedparticle.h"
#include "meshops.h"
#include "energy.h"


void ParticleAdhesion::find_pairs(Eigen::MatrixXd V, Eigen::MatrixXd V_particle, double distance_threshold, std::vector<std::pair<int, int>>& bonds, std::vector<std::pair<int, int>>& permanent_bonds)
{
    nearest.resize(V.rows());
    nearest.setConstant(-1);

    // Copy the existing permanent bonds to the new bonds list
    bonds = permanent_bonds;

    for (int i = 0; i < V.rows(); i++) {
        double min_distance = std::numeric_limits<double>::max();
        int nearest_index = -1;

        for (int j = 0; j < V_particle.rows(); j++) {
            double distance = sqrt(pow(V(i, 0) - V_particle(j, 0), 2)
                + pow(V(i, 1) - V_particle(j, 1), 2)
                + pow(V(i, 2) - V_particle(j, 2), 2));

            if (distance < min_distance) {
                min_distance = distance;
                nearest_index = j;
            }
        }
      
        if (nearest_index != -1 && min_distance < distance_threshold) {
            // Check if the bond (i, nearest_index) already exists in bonds or permanent_bonds
            bool bond_exists = false;
            for (const auto& bond : bonds) {
                if ((bond.first == i && bond.second == nearest_index) ||
                    (bond.first == nearest_index && bond.second == i)) {
                    bond_exists = true;
                    break;
                }
            }

            for (const auto& bond : permanent_bonds) {
                if ((bond.first == i && bond.second == nearest_index) ||
                    (bond.first == nearest_index && bond.second == i)) {
                    bond_exists = true;
                    break;
                }
            }

            if (!bond_exists) {
                // If the bond doesn't exist in either list, add it to the bonds vector
                bonds.emplace_back(i, nearest_index);
            }
        }
    }

    // Make the current bonds list permanent for the next iteration
    permanent_bonds = bonds;

std::ofstream outputFile("bonds.txt");

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    for (const auto& bond : bonds) {
        int vertex1_idx = bond.first;
        int vertex2_idx = bond.second;

        // Calculate the bond length between the two vertices
        double distance = (V.row(vertex1_idx) - V_particle.row(vertex2_idx)).norm();

        // Write the bond information to the file
        outputFile << "Bond: (" << vertex1_idx << ", " << vertex2_idx << "), Distance: " << distance << std::endl;
    }

    outputFile.close();
}


// void ParticleAdhesion::find_pairs(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd V_particle, Eigen::MatrixXi F_particle,double distance_threshold, std::vector<std::pair<int, int>> &bonds)
// {
// //Mesh M1;
//  //const double distance_threshold = 0.20;
// nearest.resize(V.rows());
// nearest.setConstant(-1); // initialize to -1 to mark vertices that don't have a nearest neighbor
// for (int i = 0; i < V.rows(); i++) {
//     double min_distance = std::numeric_limits<double>::max();
//     int nearest_index = -1;
//     for (int j = 0; j < V_particle.rows(); j++) {
//         double distance = sqrt(pow(V(i,0)-V_particle(j,0), 2)
//                              +pow(V(i,1)-V_particle(j,1), 2)
//                              +pow(V(i,2)-V_particle(j,2), 2));
//         if (distance < min_distance) {
//             min_distance = distance;
//             nearest_index = j;
//         }
//     }
//     if (nearest_index != -1 && min_distance < distance_threshold) {
//         nearest(i) = nearest_index;
//     }
// }
// //std::cout << nearest<< std::endl;
// // Create bond connections based on nearest neighbors
// //std::vector<std::pair<int, int>> bonds;
// for (int i = 0; i < V.rows(); i++) {
//     int j = nearest(i);
//     if (j != -1) {
//         bonds.emplace_back(i, j);
//     }
// }

// Output the number of bonds and the list of bonds
// std::cout << "Number of bonds: " << bonds.size() << std::endl;
// for (const auto& bond : bonds) {
//     std::cout << bond.first << " " << bond.second << std::endl;
// }   
//return bonds;
//}


void ParticleAdhesion::remove_long_bonds(std::vector<std::pair<int, int>>& bonds, Eigen::MatrixXd& V, Eigen::MatrixXd& V_particle, double max_bond_length)
{
    // Create a new vector to store the updated bonds
    std::vector<std::pair<int, int>> updated_bonds;

    for (const auto& bond : bonds) {
        int vertex1_idx = bond.first;
        int vertex2_idx = bond.second;

        // Calculate the bond length between the two vertices
        double bond_length = (V.row(vertex1_idx) - V_particle.row(vertex2_idx)).norm();

        // Check if the bond length is less than or equal to the maximum allowed length
        if (bond_length <= max_bond_length) {
            // If the bond length is within the threshold, keep this bond
            updated_bonds.push_back(bond);
        }
        // If the bond length is greater than the threshold, skip this bond (delete it).
    }

    // Replace the original bonds vector with the updated one
    bonds = updated_bonds;
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