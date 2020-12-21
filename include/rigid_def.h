//
// Created by behrooz on 12/14/20.
//

#ifndef RIGID_DEF_H
#define RIGID_DEF_H
#include <igl/unproject.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
//std
#include <vector>
#include <string>

//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

#include <iostream>
#include <vector>

#include <visualization.h>

namespace rigid_body{
    class Boundary{
    public:
        Boundary() = default;
        ~Boundary() = default;
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
        void init_boundary_2d(Eigen::MatrixXd &q,
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

            // Add edges and add it to the geometry for further boundary checking
            int obst_id = Visualize::viewer().append_mesh();
            geometry_id.push_back(std::pair<int, std::string> (obst_id, "boundary"));
            geometry.push_back(std::pair<Eigen::MatrixXd, Eigen::MatrixXi> (box_V, box_F));

            for (unsigned i=0; i < box_F.rows(); i++){
                Visualize::viewer().data_list[obst_id].add_edges(
                        box_V.row(box_F(i,0)),
                        box_V.row(box_F(i,1)),
                        Eigen::RowVector3d(1,0,0)
                );
            }
            Visualize::viewer().core().camera_zoom = 0.1;
            Visualize::viewer().core().object_scale = 1.0;
            std::cout << particles_viewer_id << std::endl;
            std::cout << "Setup Done" << std::endl;
        }

    };
}



#endif //RIGID_DEF_H
