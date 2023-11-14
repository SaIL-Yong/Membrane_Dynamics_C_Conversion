#pragma once
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#endif
#include <iostream>
#include <fstream>
#include <string>

void readParameter();
struct Parameter {
    int iterations, tolerance_flag, forced_wrapping_flag;
    int logfrequency, dumpfrequency, resfrequency;
    double tolfrequency;
    int particle_flag, particle_position;
    double dt, Kb, Kv, Ka, gamma, gammap, tolerance, wrapping_fraction, wrapping_bias_strength;
    double reduced_volume, particle_radius, adhesion_strength, potential_range;
    std::string meshFile, outFile, resFile;
    int mesh_reg_frequency, vertex_smoothing_flag, delaunay_triangulation_flag;
    int angle_condition_flag;
    double X0, Y0, Z0;
    double X1,Y1,Z1;
    double C0;///spontaneous_curvature_input
    int particle_coord_flag;
};
