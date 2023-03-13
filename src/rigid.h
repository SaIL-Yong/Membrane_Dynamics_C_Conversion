#pragma once
#ifndef RIGID_H_
#define RIGID_H_
#endif
#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
#include "meshedparticle.h"
#include <igl/point_mesh_squared_distance.h>
#include <vector>
#include "meshops.h"

struct RigidBody{
    //constant quantities which remains same through out the simulation
    Eigen::Matrix3f Inertia,Invertia_inv;
    //State Variables
    Eigen::Vector3f rigid_body_center,linear_momentum,angular_momentum;
    Eigen::Quaterniond q();  
    //auxiliary variables
    Eigen::Matrix3f rotation_matrix,Iinv;
    Eigen::Vector3f velocity, omega,
    //pre-computed propertioes
    Eigen::Vector3f body_force, body_torque;

    
} ;
class RigidBodyDynamics{
public:
 void compute_body_center(Eigen::MatrixXd V, RigidBody* rigid_body) {
        // implementation to compute body center based on input matrix V
        // store the result in rigid_body->rigid_body_center
        rigid_body->rigid_body_center = Eigen::Vector3f::Zero(); // example value
 }

void compute_body_force(Eigen::MatrixXd V,Eigen::MatrixXd Force_Adhesion, RigidBody* rigid_body) {
        rigid_body->body_foce = Eigen::Vector3f::Zero(); // calculate from negative of Force_Adhesion
        
 }

void compute_body_torque(Eigen::MatrixXd V,Eigen::MatrixXd Force_Adhesion, RigidBody* rigid_body) {
        
        rigid_body->body_torque = Eigen::Vector3f::Zero(); // calculate from negative of Force_Adhesion
        
 }

void compute_angular_momentum( RigidBody* rigid_body){

}
void compute_Inetria(Eigen::MatrixXd V, RigidBody* rigid_body) {
        rigid_body->Inertia = Eigen::Matrix3f::Zero(); 
        rigid_body->Inertia_inv = Eigen::Matrix3f::Zero();
 }

void compute_rotation_matrix( RigidBody* rigid_body) {
        
        rigid_body->rotation_matrix =  q.toRotationMatrix();;
 }


void Compute_Force_and_Torque(double timestep, RigidBody *rigid_body){

}

void rigid_body_euler_simulation();



};



