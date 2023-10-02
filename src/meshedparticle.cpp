#include "meshedparticle.h"
#include "meshops.h"
#include "energy.h"


//Permanent Bond Creation
/*
void ParticleAdhesion::find_pairs(Eigen::MatrixXd V, Eigen::MatrixXd V_particle, double distance_threshold, std::vector<std::pair<int, int>>& bonds)
{
    nearest.resize(V.rows());
    nearest.setConstant(-1);

    // Create a new vector to store the updated bonds
    std::vector<std::pair<int, int>> updated_bonds = bonds;

    for (int i = 0; i < V.rows(); i++) {
        bool membv_occupied = false;
        for (const auto& bond : updated_bonds) {
            if (bond.first == i) {
                membv_occupied = true;
                break;
            }
        }
        if (membv_occupied) continue;

        double min_distance = std::numeric_limits<double>::max();
        int nearest_index = -1;

        for (int j = 0; j < V_particle.rows(); j++) {
            bool particlev_occupied = false;
            for (const auto& bond : updated_bonds) {
                if (bond.second == j) {
                particlev_occupied = true;
                break;
                }
            }
            if (particlev_occupied) continue;
            
            double distance = sqrt(pow(V(i, 0) - V_particle(j, 0), 2)
                + pow(V(i, 1) - V_particle(j, 1), 2)
                + pow(V(i, 2) - V_particle(j, 2), 2));

            if (distance < min_distance) {
                min_distance = distance;
                nearest_index = j;
            }
        }

        if (nearest_index != -1 && min_distance < distance_threshold) updated_bonds.emplace_back(i, nearest_index);

    }

    // Update the bonds vector with the updated_bonds
    bonds = updated_bonds;


    std::ofstream outputFile("bonds.txt");

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    for (const auto& bond : bonds) {
        int vertex1_idx = bond.first;
        int vertex2_idx = bond.second;

         // Write the bond information to the file
        outputFile << "Bond: (" << vertex1_idx << ", " << vertex2_idx << "), Distance: " << distance << std::endl;
    }

    outputFile.close();
}

*/


//Temporary Bond Creation
void ParticleAdhesion::find_pairs(Eigen::MatrixXd V,Eigen::MatrixXi F, Eigen::MatrixXd V_particle, double distance_threshold, std::vector<std::pair<int, int>>& bonds)
{
    bonds.clear(); // Clear any existing bonds since you don't want to keep them

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
            bonds.emplace_back(i, nearest_index); // Generate a new bond based on distance
        }
    }

    std::ofstream outputFile("bonds.txt");

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    for (const auto& bond : bonds) {
        int vertex1_idx = bond.first;
        int vertex2_idx = bond.second;

        // Write the bond information to the file
        double distance = sqrt(pow(V(vertex1_idx, 0) - V_particle(vertex2_idx, 0), 2)
            + pow(V(vertex1_idx, 1) - V_particle(vertex2_idx, 1), 2)
            + pow(V(vertex1_idx, 2) - V_particle(vertex2_idx, 2), 2));

        outputFile << "Bond: (" << vertex1_idx << ", " << vertex2_idx << "), Distance: " << distance << std::endl;
    }

    outputFile.close();
    // Eigen::MatrixXd P = Eigen::MatrixXd::Random(100,3);
    // Eigen::VectorXd sqrD;
    // Eigen::VectorXi I;
    // Eigen::MatrixXd C;
    // igl::point_mesh_squared_distance(V_particle,V,F,sqrD,I,C);
    

    // std::cout << "I: " << I << std::endl;
    // std::cout << "C: " << C << std::endl;
    // std::cout << "sqrD: " << sqrD << std::endl;

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

