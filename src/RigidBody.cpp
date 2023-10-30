//  #include "RigidBody.h"
// include <Eigen/Core>
// #include <igl/massmatrix.h>
// #include <igl/doublearea.h>


// // Convert the state of a rigid body with a triangulated mesh to an array
// void State_to_Array(const RigidBody& rb, double* y) {
//     // Copy position components
//     *y++ = rb.x[0];
//     *y++ = rb.x[1];
//     *y++ = rb.x[2];

//     // Copy rotation matrix components
//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             *y++ = rb.R(i, j);
//         }
//     }

//     // Copy linear momentum components
//     *y++ = rb.P[0];
//     *y++ = rb.P[1];
//     *y++ = rb.P[2];

//     // Copy angular momentum components
//     *y++ = rb.L[0];
//     *y++ = rb.L[1];
//     *y++ = rb.L[2];

//     // You may need to add code here to copy the mesh information (vertices and faces) if needed.
// }

// // Convert an array to the state of a rigid body with a triangulated mesh
// void Array_to_State(RigidBody& rb, const double* y) {
//     // Copy position components
//     rb.x[0] = *y++;
//     rb.x[1] = *y++;
//     rb.x[2] = *y++;

//     // Copy rotation matrix components
//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             rb.R(i, j) = *y++;
//         }
//     }

//     // Copy linear momentum components
//     rb.P[0] = *y++;
//     rb.P[1] = *y++;
//     rb.P[2] = *y++;

//     // Copy angular momentum components
//     rb.L[0] = *y++;
//     rb.L[1] = *y++;
//     rb.L[2] = *y++;

//     // You may need to add code here to copy the mesh information (vertices and faces) if needed.

//     // Recalculate auxiliary variables
//     rb.v = rb.P / rb.mass;
//     rb.Iinv = rb.R * rb.Ibodyinv * rb.R.transpose();
//     rb.omega = rb.Iinv * rb.L;
// }

// void ddt_State_to_Array(RigidBody *rb, double *ydot) {
//     /* Copy v(t) = ẋ(t) into ydot */
//     *ydot++ = rb->v[0];
//     *ydot++ = rb->v[1];
//     *ydot++ = rb->v[2];

//     /* Compute Ṙ(t) = ω(t) * R(t) */
//     Eigen::Matrix3d Rdot = Star(rb->omega) * rb->R;

//     /* Copy Ṙ(t) into array */
//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             *ydot++ = Rdot(i, j);
//         }
//     }

//     /* Copy force components into ydot */
//     *ydot++ = rb->force[0];
//     *ydot++ = rb->force[1];
//     *ydot++ = rb->force[2];

//     /* Copy torque components into ydot */
//     *ydot++ = rb->torque[0];
//     *ydot++ = rb->torque[1];
//     *ydot++ = rb->torque[2];
// }




// // using namespace Eigen;

// // MatrixXd V; // Vertex positions
// // MatrixXi F; // Triangle indices

// // double density = 1.0; // Density of the material
// // double m; // Mass of the object
// // Vector3d center_of_mass; // Center of mass of the object
// // Matrix3d I; // Inertia tensor of the object

// // // Read in the mesh
// // igl::readOBJ("mesh.obj", V, F);

// // // Compute the mass properties
// // double volume = igl::volume(V, F) / 3.0; // Volume of the object
// // m = density * volume; // Compute the mass of the object
// // igl::centroid(V, F, center_of_mass); // Compute the center of mass of the object

// // // Compute the inertia tensor
// // I.setZero(); // Initialize the inertia tensor to zero
// // for (int i = 0; i < F.rows(); i++) {
// //     Vector3d v0 = V.row(F(i, 0));
// //     Vector3d v1 = V.row(F(i, 1));
// //     Vector3d v2 = V.row(F(i, 2));
// //     Vector3d c = (v0 + v1 + v2) / 3.0 - center_of_mass;
// //     double area = igl::doublearea(v0, v1, v2);
// //     I(0,0) += area * (c(1)*c(1) + c(2)*c(2));
// //     I(1,1) += area * (c(0)*c(0) + c(2)*c(2));
// //     I(2,2) += area * (c(0)*c(0) + c(1)*c(1));
// //     I(0,1) -= area * c(0) * c(1);
// //     I(0,2) -= area * c(0) * c(2);
// //     I(1,2) -= area * c(1) * c(2);
// // }
// // I(1,0) = I(0,1);
// // I(2,0) = I(0,2);
// // I(2,1) = I(1,2);
// // I *= density / 60.0;

// // // Output the inertia tensor
// // std::cout << "Inertia tensor:" << std::endl << I << std::endl;

// // Eigen::Quaterniond RigidBody::q()
// // {
// //     return Eigen::Quaterniond();
// // }
