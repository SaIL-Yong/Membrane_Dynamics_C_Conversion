#include <iostream>
#include <Eigen/Core>
#include "energy.h"
#include "meshedparticle.h"
#include "parameters.h"
Parameter parameter;
int main()
{
    // Load the two sets of vertices
    readParameter();
    Eigen::MatrixXd V1, V2; // Input vertices
    Eigen::MatrixXi F1, F2; // Input faces (not used in this example)
    int numF = F1.rows();
    int numV = V1.rows();
    Eigen::VectorXi nearest;                // Nearest neighbor of each vertex in V1 in V2
    std::vector<std::pair<int, int>> bonds; // Bond connections between V1 and V2

    // Load input meshes from OFF files
    igl::readOFF("equilibrated.off", V1, F1);
    igl::readOFF("particle.off", V2, F2);
    Eigen::MatrixXd Force_Area(numV, 3), Force_Volume(numV, 3), Force_Bending(numV, 3), Force_Adhesion(numV, 3), velocity(numV, 3), Force_Total(numV, 3); // force components
    velocity.setZero();
    double EnergyVolume = 0.0, EnergyArea = 0.0, EnergyBending = 0.0, EnergyAdhesion = 0.0, EnergyBias = 0.0,
           EnergyTotal = 0.0, EnergyTotalold_log = 0.0, EnergyChangeRate_log = 0.0, EnergyChangeRate_avg = 0.0;
    double Kb = 0.01;
    // parameters for particle adhesion
    int particle_flag = parameter.particle_flag;
    int particle_position = parameter.particle_position;
    double Rp, u, U, rho, rc, X0, Y0, Z0, Ew_t, Kw;
    int angle_flag;

    Rp = parameter.particle_radius;
    u = parameter.adhesion_strength;
    U = (Kb * u) / (Rp * Rp);
    rho = parameter.potential_range;
    rc = 5.0 * rho;
    angle_flag = parameter.angle_condition_flag;

    X0 = 0.0, Y0 = 0.0, Z0 = V1.col(2).maxCoeff() + parameter.particle_position * (Rp + 1.0 * rho);
    
    // double Kb = parameter.Kb;
    double rVol_t = parameter.reduced_volume;
    double Kv = 0.0;
    double Ka = parameter.Ka;
    double Rv = 1.0;
    double area_target = 4 * PI * Rv * Rv;
    Kv = parameter.Kv;
    double volume_target = rVol_t * (4.0 / 3.0) * PI * pow(Rv, 3);
    

    double distance_threshold = 0.20;
    // Calculate the distances between each pair of vertices
    ParticleAdhesion P1;
    P1.find_pairs(V1,F1, V2,F2, distance_threshold, bonds);
    std::cout << V1(bonds[0].first, 1) << std::endl;
    Mesh M1;
    Energy E1;
    M1.mesh_cal(V1, F1);
    E1.compute_bendingenergy_force(V1, F1, Kb, Force_Bending, EnergyBending, M1);
    E1.compute_areaenergy_force(V1, F1, Ka, area_target, Force_Area, EnergyArea, M1);
    E1.compute_volumeenergy_force(V1, F1, Kv, volume_target, Force_Volume, EnergyVolume, M1);
    E1.compute_adhesion_energy_force(V1, F1, V2, F2, rho, U, Force_Adhesion, bonds, EnergyAdhesion, M1);
    //std::cout << "Adhesion Force:" << Force_Adhesion << std::endl;
    //std::cout << "Adhesion Energy:" << EnergyAdhesion << std::endl;
}

void readParameter()
{
    std::string line;
    std::ifstream runfile;
    runfile.open("run_parameters.txt");
    if (!runfile.is_open())
        std::cout << "ERROR: cannot access simlation parameter file." << std::endl;
    getline(runfile, line);
    runfile >> parameter.iterations;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dt;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.Kb;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.Ka;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.Kv;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.reduced_volume;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.tolerance;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.tolerance_flag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.tolfrequency; // in units of sim. time (iteration * timestep)
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.gamma;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.particle_flag;
    getline(runfile, line);
    if (parameter.particle_flag)
    {
        getline(runfile, line);
        runfile >> parameter.particle_position;
        getline(runfile, line);
        getline(runfile, line);
        if (line.compare("particle_coordinate") == 0)
        {
            runfile >> parameter.X0 >> parameter.Y0 >> parameter.Z0;
            parameter.particle_coord_flag = 1;
            getline(runfile, line);
            getline(runfile, line);
        }
        else
            parameter.particle_coord_flag = 0;
        runfile >> parameter.particle_radius;
        getline(runfile, line);
        getline(runfile, line);
        runfile >> parameter.adhesion_strength;
        getline(runfile, line);
        getline(runfile, line);
        runfile >> parameter.potential_range;
        getline(runfile, line);
        getline(runfile, line);
        runfile >> parameter.angle_condition_flag;
        getline(runfile, line);
        getline(runfile, line);
        runfile >> parameter.forced_wrapping_flag;
        getline(runfile, line);
        getline(runfile, line);
        runfile >> parameter.wrapping_fraction;
        getline(runfile, line);
        getline(runfile, line);
        runfile >> parameter.wrapping_bias_strength;
        getline(runfile, line);
    }
    getline(runfile, line);
    getline(runfile, parameter.meshFile);
    getline(runfile, line);
    getline(runfile, parameter.particleFile);
    getline(runfile, line);
    getline(runfile, parameter.outFile);
    getline(runfile, line);
    getline(runfile, parameter.resFile);
    getline(runfile, line);
    runfile >> parameter.logfrequency;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dumpfrequency;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.resfrequency;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.mesh_reg_frequency;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.vertex_smoothing_flag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.delaunay_triangulation_flag;
    runfile.close();
}
