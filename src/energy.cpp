#include <iostream>
#include <cmath>
#include "energy.h"

void Energy::compute_bendingenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kb, Eigen::MatrixXd& Force_Bending, double& bending_energy, Mesh m)
{
  bending_energy = 0.0;
  Force_Bending.setZero();
  EB = 2.0 * Kb * (m.H_squared.transpose() * m.area_voronoi).diagonal();
  bending_energy = EB.sum();

  Lap_H = m.Minv * (m.L * m.H_signed);
  force_density = (2.0 * m.H_signed.array() * (m.H_squared - m.K).array()) + Lap_H.array();
  vector_term = force_density.array() * m.area_voronoi.array();
  Force_Bending = (2.0 *Kb)*(m.V_normals.array().colwise() * vector_term.array());
}


void Energy::compute_areaenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Ka, double area_target, Eigen::MatrixXd& Force_Area, double& area_energy, Mesh m)
{
  area_energy = 0.0;
  Force_Area.setZero();
  da = m.area_total - area_target;
  area_energy = Ka * da * da / area_target;
  scalar_term = -2.0 * Ka * da / area_target;
  AG = m.area_grad(V, F);
  Force_Area = scalar_term * AG;
}


void Energy::compute_volumeenergy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double Kv, double volume_target, Eigen::MatrixXd& Force_Volume, double& volume_energy, Mesh m)
{
  volume_energy = 0.0;
  Force_Volume.setZero();
  if (std::abs(Kv) > EPS) {
    VG = m.volume_grad(V, F);
    dv = m.volume_total - volume_target;
    volume_energy = Kv * dv * dv / volume_target;
    scalar_term = -2.0 * Kv * dv / volume_target;
    Force_Volume = scalar_term * VG;
  }
}
void Energy::compute_adhesion_energy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd V_particle, Eigen::MatrixXi F_particle,double rho, double U,
                                           Eigen::MatrixXd& Force_Adhesion,std::vector<std::pair<int, int>> bonds, double& EnergyAdhesion,  Mesh m)
{
  Force_Adhesion.setZero();
  coefficient.resize(V.rows());
  coefficient_derivative_x.resize(V.rows());
  coefficient_derivative_y.resize(V.rows());
  coefficient_derivative_z.resize(V.rows());
  distance.resize(V.rows());
  dc.resize(V.rows());
  //Mod_Bias.resize(V.rows());
  coefficient.setZero();
  coefficient_derivative_x.setZero();
  coefficient_derivative_y.setZero();
  coefficient_derivative_z.setZero();
  distance.setZero();
  dc.setZero();


  // vector connecting cetner of the particle and to the vertices;
  //comvec = V.rowwise() - particle_center;
  double distance_threshold = 0.20;
  //std::vector<std::pair<int, int>> bonds; 
  ParticleAdhesion P1;
  P1.find_pairs(V,F,V_particle,F_particle,distance_threshold,bonds);

  // can this loop be replaced by the eigen operation?
  for (int i = 0; i <bonds.size(); i++) {
    //distance(i) = sqrt((V(i,0)-X)*(V(i,0)-X)+(V(i,1)-Y)*(V(i,1)-Y)+(V(i,2)-Z)*(V(i,2)-Z));
    dc(bonds[i].first) = sqrt((V(bonds[i].first,0)-(V_particle(bonds[i].second,0)))*(V(bonds[i].first,0)-(V_particle(bonds[i].second,0)))+
    (V(bonds[i].first,1)-(V_particle(bonds[i].second,1)))*(V(bonds[i].first,1)-(V_particle(bonds[i].second,1)))+
    (V(bonds[i].first,2)-(V_particle(bonds[i].second,2)))*(V(bonds[i].first,2)-(V_particle(bonds[i].second,2))));

    // angle between connecting vector and vertex normal
    //angle = acos((comvec.row(i)).dot(m.V_normals.row(i)) / (comvec.row(i).norm() * m.V_normals.row(i).norm()));
    
    // if (std::abs(dc(i)) > rc) continue;
    // if (angle_flag) {
    //   if (particle_position > 0 && angle <= 0.5*PI) continue;
    //   if (particle_position < 0 && angle >= 0.5*PI) continue;
    // }

    coefficient.row(bonds[i].first) << U * (exp(-(2.0*dc(bonds[i].first))/rho) - 2.0*(exp(-dc(bonds[i].first)/rho)));
    coefficient_derivative_x.row(bonds[i].first) << (U/dc(bonds[i].first)*rho)
                                *(-exp(-(2.0*dc(bonds[i].first))/rho) + exp(-dc(bonds[i].first)/rho)) * 2.0 * (V(bonds[i].first,0)-(V_particle(bonds[i].second,0)));
    coefficient_derivative_y.row(bonds[i].first) << (U/(dc(bonds[i].first)*rho))
                                *(-exp(-(2.0*dc(bonds[i].first))/rho) + exp(-dc(bonds[i].first)/rho)) * 2.0 * (V(bonds[i].first,1)-(V_particle(bonds[i].second,1)));
    coefficient_derivative_z.row(bonds[i].first) << (U/dc(bonds[i].first)*rho)
                                *(-exp(-(2.0*dc(bonds[i].first))/rho) + exp(-dc(bonds[i].first)/rho)) * 2.0 * (V(bonds[i].first,2)-(V_particle(bonds[i].second,2)));

    // if (dc(i) > EPS && std::abs(Kw) > EPS) Mod_Bias(i) = 1.0;
  }
  //std::cout<<coefficient<<std::endl;
  coefficient_of_derivative.resize(V.rows(), 3);
  coefficient_of_derivative.col(0)=coefficient_derivative_x.transpose();
  coefficient_of_derivative.col(1)=coefficient_derivative_y.transpose();
  coefficient_of_derivative.col(2)=coefficient_derivative_z.transpose();
  
  First_Term = -(AG.array().colwise()*coefficient.array());
  Second_Term = -(coefficient_of_derivative.array().colwise()*m.area_voronoi.array());

  Sum = First_Term + Second_Term;
  Ead = coefficient.array()*m.area_voronoi.array();

  EnergyAdhesion = Ead.sum();

  Force_Adhesion = Sum; 
}

// void Energy::compute_adhesion_energy_force(Eigen::MatrixXd V, Eigen::MatrixXi F, double X, double Y, double Z,
//                                            double Rp, double rho, double U, double rc, int angle_flag, int particle_position, double Ew_t, double Kw,
//                                            Eigen::MatrixXd& Force_Adhesion, double& EnergyAdhesion, double& EnergyBias, Mesh m)
// {
//   Force_Adhesion.setZero();
//   coefficient.resize(V.rows());
//   coefficient_derivative_x.resize(V.rows());
//   coefficient_derivative_y.resize(V.rows());
//   coefficient_derivative_z.resize(V.rows());
//   distance.resize(V.rows());
//   dc.resize(V.rows());
//   //Mod_Bias.resize(V.rows());
//   coefficient.setZero();
//   coefficient_derivative_x.setZero();
//   coefficient_derivative_y.setZero();
//   coefficient_derivative_z.setZero();
//   distance.setZero();
//   dc.setZero();
//   //Mod_Bias.setZero();
//   EnergyAdhesion = 0.0;
//   EnergyBias = 0.0;

//   particle_center<<X, Y, Z;
//   // vector connecting cetner of the particle and to the vertices;
//   comvec = V.rowwise() - particle_center;

//   // can this loop be replaced by the eigen operation?
//   for (int i = 0; i < V.rows(); i++) {
//     distance(i) = sqrt((V(i,0)-X)*(V(i,0)-X)+(V(i,1)-Y)*(V(i,1)-Y)+(V(i,2)-Z)*(V(i,2)-Z));
//     dc(i) = distance(i) - Rp;

//     // angle between connecting vector and vertex normal
//     angle = acos((comvec.row(i)).dot(m.V_normals.row(i)) / (comvec.row(i).norm() * m.V_normals.row(i).norm()));
    
//     if (std::abs(dc(i)) > rc) continue;
//     if (angle_flag) {
//       if (particle_position > 0 && angle <= 0.5*PI) continue;
//       if (particle_position < 0 && angle >= 0.5*PI) continue;
//     }

//     coefficient(i) = U * (exp(-(2.0*dc(i))/rho) - 2.0*exp(-dc(i)/rho));
//     coefficient_derivative_x(i) = (U/(distance(i)*rho))
//                                 *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,0)-X);
//     coefficient_derivative_y(i) = (U/(distance(i)*rho))
//                                 *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,1)-Y);
//     coefficient_derivative_z(i) = (U/(distance(i)*rho))
//                                 *(-exp(-(2.0*dc(i))/rho) + exp(-dc(i)/rho)) * 2.0 * (V(i,2)-Z);

//     // if (dc(i) > EPS && std::abs(Kw) > EPS) Mod_Bias(i) = 1.0;
//   }

//   coefficient_of_derivative.resize(V.rows(), 3);
//   coefficient_of_derivative.col(0)=coefficient_derivative_x.transpose();
//   coefficient_of_derivative.col(1)=coefficient_derivative_y.transpose();
//   coefficient_of_derivative.col(2)=coefficient_derivative_z.transpose();
  
//   First_Term = -(AG.array().colwise()*coefficient.array());
//   Second_Term = -(coefficient_of_derivative.array().colwise()*m.area_voronoi.array());

//   Sum = First_Term + Second_Term;
//   Ead = coefficient.array()*m.area_voronoi.array();

//   EnergyAdhesion = Ead.sum();

//   if (std::abs(Kw) > EPS) {
//     dEw = EnergyAdhesion - Ew_t;
//     EnergyBias = 0.5 * Kw * dEw * dEw;
//     //Force_Biased = Kw * dEw * (Sum.array().colwise() * Mod_Bias.array());
//     Force_Biased = Kw * dEw * Sum;
//     Force_Adhesion = Sum + Force_Biased;
//   } else Force_Adhesion = Sum; 
// }