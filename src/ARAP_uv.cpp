// author: @janpiotraschke
// date: 2023-01-11
// license: Apache License 2.0
// description: ARAP (As-Rigid-As-Possible) parameterization
// g++ -std=c++11 -lpthread -I ./libigl/include/ -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 src/ARAP_uv.cpp -o src/ARAP_uv


// ! look at https://github.com/libigl/libigl/tree/main/include/igl to see all possible headers
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOFF.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <iostream>
#include <filesystem>

// ? bounding: https://github.com/libigl/libigl/blob/main/include/igl/bounding_box.h

Eigen::MatrixXd vertices;  // double, dynamic matrix
Eigen::MatrixXi faces;  // integer, dynamic matrix
Eigen::MatrixXd vertices_uv;
Eigen::MatrixXd initial_guess;  // list of weights


int main(int argc, char *argv[])
{
    std::__fs::filesystem::path cwd = std::__fs::filesystem::current_path();
    std::cout << cwd << std::endl;  // print current working directory
    // TODO: change to relative path
    // igl::read_triangle_mesh("git_repos/Confined_active_particles/meshes/ellipsoid_x4.stl", vertices, faces);  // Load a mesh in STL format:
    igl::readOFF("git_repos/Confined_active_particles/meshes/camelhead.off", vertices, faces);  // Load a mesh in OFF format:

    // Compute the initial solution for ARAP (harmonic parametrization)
    Eigen::VectorXi bnd;
    igl::boundary_loop(faces, bnd);
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(vertices, bnd, bnd_uv);

    igl::harmonic(vertices, faces, bnd, bnd_uv, 1, initial_guess);  // the '1' is the k power of harmonic operation (1: harmonic, 2: biharmonic, etc)

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;
    Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

    // Initialize ARAP
    arap_data.max_iter = 100;
    igl::arap_precomputation(vertices, faces, 2, b, arap_data);  // '2' means that we're going to *solve* in 2d

    auto vertices_uv = initial_guess;  // Solve arap using the harmonic map as initial guess

    igl::arap_solve(bc, arap_data, vertices_uv);  // parameterizes the mesh

    vertices_uv *= 20;  // Scale UV to make the texture more clear

    // the following two lines are for the usage of the igl::write_triangle_mesh function
    vertices_uv.conservativeResize(Eigen::NoChange, vertices_uv.cols()+1);
    vertices_uv.col(vertices_uv.cols()-1).setZero();

    igl::write_triangle_mesh("git_repos/Confined_active_particles/meshes/camel_uv.stl", vertices_uv, faces);  // write the mesh in STL format

    return 0;
}
