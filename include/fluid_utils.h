//
// Created by behrooz on 12/14/20.
//

#ifndef FLUID_UTILS_H
#define FLUID_UTILS_H

#include <visualization.h>
//std
#include <vector>
#include <string>

#include <visualization.h>
//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <fluid_utils.h>
#include <iostream>
#include <vector>

namespace fluid{
    // ========================== Particle stuff =================================
    // Useful Variables
    int x_size_max;//The x axis in meter
    int y_size_max;//The y axis in meter
    int x_size_min;//The x axis in meter
    int y_size_min;//The y axis in meter
    int x_size_range;
    int y_size_range;

    int particles_viewer_id;//The id that we can use to access particles list in the viewer to update them
    int num_particles; //number of particles inside the simulation

    /*
     * @brief This function create a rectangle and add some particles inside this rectangle for basic simulation
     * @param Viewer: is the opengl viewer of the scene
     * @param q: is the general coordination for only particles
     * @param qdot: is the velocity of particles
     * @param geometry: is the buffer for storing all the obstacles inside our simulation scene
     * @param geometry_id: is the buffer to store the Viewer.data_list ids for obstacles
     */
    void init_fluid_2d(igl::opengl::glfw::Viewer &Viewer,
                       Eigen::MatrixXd &q,
                       Eigen::MatrixXd &qdot,
                       std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> &geometry,
                       std::vector<std::pair<int, std::string>> &geometry_id){
        std::cout << "Setting up the simulation particles" << std::endl;
        // Create The boundary box of the simulation (a square)
        //Set simulation size
        int x_min = -10;
        int x_max = 10;
        int y_min = -10;
        int y_max = 10;

        int x_range = x_max - x_min;
        int y_range = y_max - y_min;

        x_size_max = x_max;
        y_size_max = y_max;
        x_size_min = x_min;
        y_size_min = y_min;
        x_size_range = x_range;
        y_size_range = y_range;

        //Vertex of bounding box
        Eigen::MatrixXd box_V(4,3);
        box_V <<
            x_min, y_min, 0,
            x_max, y_min, 0,
            x_max, y_max, 0,
            x_min, y_max, 0;
        // Edges of bounding box
        Eigen::MatrixXi box_F(4,2);
        box_F <<
            0, 1,
            1, 2,
            2, 3,
            3, 0;

        // Create The Fluid inside the box (a rectangle shape points)
        num_particles = 100;
        // Start to create the particles in the top left of the square
        double fluid_x_range = x_size_range / 8.0;
        double fluid_y_range = y_size_range / 2.0;
        // fill the rectangle with 100 particles
        // So be careful about the initial number of particles it should have int sqrt.
        int col_num_particles = sqrt(num_particles / 4);
        int row_num_particles = col_num_particles * 4;
        int offset = 1;// For not trapping in the initialize boundary conditions
        double x_space = fluid_x_range / col_num_particles;
        double y_space = fluid_y_range / row_num_particles;
        // I want it to fill the left of the rectangle.
        q.resize(num_particles,3);
        qdot.resize(num_particles, 3);
        qdot.setZero();
        int cnt = 0;
        for(int req_row = 0; req_row < row_num_particles; req_row++){
            for(int req_col = 0; req_col < col_num_particles; req_col++){
                Eigen::Vector3d pos;
                pos(0) = x_size_min + req_col * (x_space) + offset;
                pos(1) = y_size_max - req_row * (y_space) - offset;
                pos(2) = 0;
                q.row(cnt) = pos;
                qdot.row(cnt).setZero();
                cnt++;
            }
        }
        assert(num_particles == cnt);
        if(num_particles != cnt){
            std::cout << "The number of requested particles is"
                         " not equal to the number of particles in the simulation" << std::endl;
        }
        particles_viewer_id = Viewer.append_mesh();
        Viewer.data_list[particles_viewer_id].set_points(q, Eigen::RowVector3d(0,0,128));
        Viewer.data().point_size = 5;

        // Add edges and add it to the geometry for further boundary checking
        int obst_id = Viewer.append_mesh();
        geometry_id.push_back(std::pair<int, std::string> (obst_id, "boundary"));
        geometry.push_back(std::pair<Eigen::MatrixXd, Eigen::MatrixXi> (box_V, box_F));

        for (unsigned i=0; i < box_F.rows(); i++){
            Viewer.data_list[obst_id].add_edges(
                            box_V.row(box_F(i,0)),
                            box_V.row(box_F(i,1)),
                            Eigen::RowVector3d(1,0,0)
                    );
        }
        Viewer.core().camera_zoom = 0.1;
        Viewer.core().object_scale = 1.0;
        std::cout << "Setup Done" << std::endl;
    }

    int get_particels_data_list_id() {return particles_viewer_id;}
    // ========================== Grid stuff =================================
    int x_grid_size;// Number of cell in x direction
    int y_grid_size;// Number of cell in y direction
    //Test simulation - only applying gravity

}



#endif //FLUID_UTILS_H
