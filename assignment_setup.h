#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H
//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
//GUI
#include <visualization.h>
//STL
#include <string>
#include <iostream>
//simulation headers
#include <fluid_utils.h>
#include <simulation.h>
#include <grid.h>



//obstacle body geometry
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> geometry;
std::vector<std::pair<int, std::string>> geometry_id;


Eigen::Vector3d g(0, -9.8, 0);
Eigen::VectorXd forces;

bool simulation_pause = true;
//selection spring
double k_selected = 1e-1;

inline void updateFluidParticles(const Eigen::MatrixXd &q){
    int particle_idx = fluid::get_particels_data_list_id();
    Visualize::viewer().data_list[particle_idx].set_points(q, Eigen::RowVector3d(0,0,128));
    Visualize::viewer().data().point_size = 5;
    Visualize::viewer().data_list[particle_idx].dirty |= igl::opengl::MeshGL::DIRTY_POSITION;
}



inline void simulate(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, double dt, double t) {
    if(!simulation_pause) {
        //Simple simulation
       simulation::water_sim.simulate_one_step(q, qdot, forces, grid::grid2d);
    }
}


inline void draw(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, double t) {
    // Updating the fluid particles
    updateFluidParticles(q);
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

inline void assignment_setup(int argc, char **argv, Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, double dt) {

    // Init fluid particles and the boundary
    fluid::init_fluid_2d(Visualize::viewer(), q, qdot, geometry, geometry_id);
    //Init the Grid
    grid::grid2d.init_2d_grid(fluid::x_size_min,
                              fluid::x_size_max,
                              fluid::y_size_min,
                              fluid::y_size_max,
                              10, 10);
    grid::grid2d.show_grid(Visualize::viewer());
    // initialize the water simulation object in simulation namespace

    //Water config
    double density = 1000;
    double viscosity = 0;
    simulation::water_sim.init_simulation(dt, g, density, viscosity, geometry, geometry_id);

    // initialize the function for keyboard actions
    Visualize::viewer().callback_key_down = key_down_callback;
}

#endif

