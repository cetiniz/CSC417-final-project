//
// Created by behrooz on 12/16/20.
//

#ifndef INC_2D_GRID_H
#define INC_2D_GRID_H

#include <fluid_utils.h>
#include <visualization.h>

namespace grid{

    class Grid{
    public:
        Grid() = default;
        ~Grid() = default;
    };

    class GridDense2D: Grid{
    public:
        GridDense2D() = default;
        ~GridDense2D() = default;

        /*
         * @brief initialize the grids inside the containers
         * @param container_x_size the container size in x direction
         * @param container_y_size the container size in y direction
         * @param number of grid cell in x direction
         * @param number of grid cell in y direction
         */
        void init_2d_grid(double container_x_min,
                          double container_x_max,
                          double container_y_min,
                          double container_y_max,
                          int num_cell_x, int num_cell_y){
            this-> container_x_max = container_x_max;
            this-> container_x_min = container_x_min;
            this-> container_y_max = container_y_max;
            this-> container_y_min = container_y_min;
            this-> Nx = num_cell_x;
            this-> Ny = num_cell_y;
            this-> N = Nx * Ny;
            this-> pressure.resize(Nx * Ny);
            this-> pressure.setZero();
            this-> ux.resize((Nx + 1) * Ny);
            this-> ux.setZero();
            this-> uy.resize(Nx * (Ny + 1));
            this-> uy.setZero();
            this->cell_size_x = (container_x_max - container_x_min) / Nx;
            this->cell_size_y = (container_y_max - container_y_min) / Ny;
            //Defining Flags
            this->fluid_cell.resize(Nx * Ny);
            this->air_cell.resize(Nx * Ny);
            this->solid_cell.resize(Nx * Ny);
        }

        /*
         * Show the grid inside the viewer window
         */
        void show_grid(igl::opengl::glfw::Viewer &Viewer){
            // Horizontal grid lines
            for (unsigned y=1; y < Ny; y++){
                Viewer.data().add_edges(
                        Eigen::RowVector3d(container_x_min,
                                           y * cell_size_y + container_y_min,
                                           0),
                        Eigen::RowVector3d(container_x_max,
                                           y * cell_size_y + container_y_min,
                                           0),
                        Eigen::RowVector3d(0,255,0)
                );
            }
            // Vertical grid lines
            for (unsigned x=1; x < Nx; x++){
                Viewer.data().add_edges(
                        Eigen::RowVector3d(x * cell_size_x + container_x_min,
                                           container_y_min,
                                           0),
                        Eigen::RowVector3d(x * cell_size_x + container_x_min,
                                           container_y_max,
                                           0),
                        Eigen::RowVector3d(0,255,0)
                );
            }
        }
        // Functions to work with MAC variables

    private:
        //Number of grids in each direction
        int Nx;
        int Ny;
        int N;
        //Size of grids in each direction
        double cell_size_x;
        double cell_size_y;

        //Container Coordination
        double container_x_min;
        double container_x_max;
        double container_y_min;
        double container_y_max;

        //Flags
        std::vector<bool> fluid_cell;
        std::vector<bool> air_cell;
        std::vector<bool> solid_cell;
        //The variables of MAC grids
        Eigen::VectorXd pressure;
        Eigen::VectorXd ux;
        Eigen::VectorXd uy;
        Eigen::VectorXd uv;
    };

    GridDense2D grid2d;
}
#endif //INC_2D_GRID_H
