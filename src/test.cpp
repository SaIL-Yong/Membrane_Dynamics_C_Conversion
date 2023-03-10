#include <iostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/knn.h>
#include <Eigen/Core>
#include <vector>

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::MatrixXi MatrixXi;
typedef Eigen::VectorXi VectorXi;
typedef Eigen::VectorXd VectorXd;

// Find nearest neighbor of each vertex in V1 in V2
void find_nearest_neighbors(const MatrixXd& V1, const MatrixXd& V2,
                             VectorXi& nearest) {
    igl::knn(V1, V2, 1, nearest);
}

int main()
{
    MatrixXd V1, V2; // Input vertices
    MatrixXi F1, F2; // Input faces (not used in this example)
    VectorXi nearest; // Nearest neighbor of each vertex in V1 in V2
    std::vector<std::pair<int, int>> bonds; // Bond connections between V1 and V2
    
    // Load input meshes from OFF files
    igl::readOFF("equilibrated.off", V1, F1);
    igl::readOFF("particle.off", V2, F2);
    
    // Find nearest neighbor of each vertex in V1 in V2
    find_nearest_neighbors(V1, V2, nearest);
    
    // Create bond connections between V1 and V2 based on distance
    for (int i = 0; i < V1.rows(); i++) {
        int j = nearest(i);
        double dist = (V1.row(i) - V2.row(j)).norm();
        if (dist < threshold) { // Only create bond if distance is less than threshold
            bonds.push_back(std::make_pair(i, j));
        }
    }
    
    // Save output mesh with bond connections
    igl::writeOFF("output.off", V1, F1, bonds);
    
    return 0;
}

// #include <iostream>
// #include <fstream>
// #include <chrono>
// #include <cmath>
// #include "meshops.h"
// #include "energy.h"
// #include "parameters.h"
// #include "meshedparticle.h"
// using namespace std::chrono;
// int numV;                                               // number of vertices
// int numF;                                               // number of faces
// Eigen::MatrixXd V1,V2;                                      // matrix storing vertice coordinates
// Eigen::MatrixXi F1,F2;
// Parameter parameter;

// int main()
// {   
//     readParameter();
//     // Load two meshes
//     igl::readOFF(parameter.meshFile, V1, F1);
//     igl::readOFF(parameter.particleFile, V2, F2);
//     // Find vertex correspondence based on distance
//     ParticleAdhesion P1;
//     std::vector<std::pair<int, int>> pairs;
//     double distance_threshold = 0.1;
//     P1.find_pairs(V1, V2, F1,F2,distance_threshold, pairs);
//     // Write output
//     std::ofstream out("correspondence_pairs.txt");
//     for (const auto& p : pairs) {
//         out << p.first << "" << p.second << std::endl;
//     }
//     out.close();
//     return 0;
// }
// void readParameter()
// {
//   std::string line;
//   std::ifstream runfile;
//   runfile.open("run_parameters.txt");
//   if (!runfile.is_open()) std::cout<<"ERROR: cannot access simlation parameter file."<<std::endl;
//   getline(runfile, line);
//   runfile >> parameter.iterations;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.dt;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.Kb;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.Ka;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.Kv;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.reduced_volume;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.tolerance;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.tolerance_flag;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.tolfrequency;  // in units of sim. time (iteration * timestep)
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.gamma;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.particle_flag;
//   getline(runfile, line);
//   if (parameter.particle_flag) {
//     getline(runfile, line);
//     runfile >> parameter.particle_position;
//     getline(runfile, line);
//     getline(runfile, line);
//     if (line.compare("particle_coordinate") == 0) {
//       runfile >> parameter.X0 >> parameter.Y0 >> parameter.Z0;
//       parameter.particle_coord_flag = 1;
//       getline(runfile, line);
//       getline(runfile, line);
//     }
//     else parameter.particle_coord_flag = 0;
//     runfile >> parameter.particle_radius;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.adhesion_strength;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.potential_range;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.angle_condition_flag;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.forced_wrapping_flag;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.wrapping_fraction;
//     getline(runfile, line);
//     getline(runfile, line);
//     runfile >> parameter.wrapping_bias_strength;
//     getline(runfile, line);
//   }
//   getline(runfile, line);
//   getline(runfile, parameter.meshFile);
//   getline(runfile, line);
//   getline(runfile, parameter.particleFile);
//   getline(runfile, line);
//   getline(runfile, parameter.outFile);
//   getline(runfile, line);
//   getline(runfile, parameter.resFile);
//   getline(runfile, line);
//   runfile >> parameter.logfrequency;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.dumpfrequency;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.resfrequency;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.mesh_reg_frequency;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.vertex_smoothing_flag;
//   getline(runfile, line);
//   getline(runfile, line);
//   runfile >> parameter.delaunay_triangulation_flag;
//   runfile.close();
// }
