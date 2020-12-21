#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H
//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
//STL
#include <visualization.h>
#include <string>
#include <iostream>
//simulation headers
#include <rigid_def.h>
#include <simulation.h>
#include <grid.h>


//obstacle body geometry
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> geometry;
std::vector<std::pair<int, std::string>> geometry_id;

//Water simulation object
simulation::GeneralSimulation water_sim;
rigid_body::Boundary boundary;

Eigen::Vector3d g(0, 9.8, 0);
Eigen::VectorXd forces;

bool simulation_pause = true;
//selection spring
double k_selected = 1e-1;


inline void simulate(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, double dt, double t) {
    water_sim.simulate_one_step(q, qdot, forces);
    if(!simulation_pause) {
        //Simple simulation
//        std::cout << "WEIRD" << std::endl;
//        water_sim.simulate_one_step(q, qdot, forces);
    }
}


inline void draw(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, double t) {
    //Update the particles
    Visualize::viewer().data_list[water_sim.get_prtcl_veiwer_id()].set_points(q, Eigen::RowVector3d(0,0,128));
    Visualize::viewer().data().point_size = 5;
    //Update grid types
    water_sim.get_grid().show_cell_type();
}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int moddifiers) {
    if(key == 'R') {
        //reset the simulation
    }

    if(key == 'S') {
        simulation_pause = !simulation_pause;
    }

    return false;
}

inline void assignment_setup(int argc, char **argv, Eigen::MatrixXd &q,
                             Eigen::MatrixXd &qdot, double dt) {

    // Init fluid particles and the boundary
    boundary.init_boundary_2d(q, qdot, geometry, geometry_id);

    // initialize the water simulation object in simulation namespace

    // For each rigid body in the scene
    Eigen::SparseMatrixd N;
    for (auto& geo : geometry) {
        N.resize(geo.first.rows(), geo.first.rows());
        N.setIdentity();
        Visualize::add_object_to_scene(
                geo.first,
                geo.second,
                geo.first,
                geo.second,
                N,
                Eigen::RowVector3d(244, 165, 130) / 255.
        );
    }

    //Water config
    double density = 1000;
    double viscosity = 0;
    water_sim.init_simulation(dt, g, density, viscosity, geometry, geometry_id, q, qdot);
    water_sim.get_grid().show_grid();
    // initialize the function for keyboard actions
    Visualize::viewer().callback_key_down = key_down_callback;
}

#endif

