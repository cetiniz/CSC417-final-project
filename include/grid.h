//
// Created by behrooz on 12/16/20.
//

#ifndef INC_2D_GRID_H
#define INC_2D_GRID_H

#include <rigid_def.h>
#include <visualization.h>
#include <fstream>
#include <string>
#include <math.h>

namespace grid{
    /*
      * Different types for the cells inside a macgrid
    */
    enum CellTypes {
        SOLID,
        FLUID,
        AIR
    };

    struct Cell {
        CellTypes type;
        int layer;
        std::vector<Eigen::Vector3i> neighbours;
        double pressure;
        Eigen::Vector3d v_;
        Eigen::Vector3d tmpV_; // For flip
        Eigen::Vector3i coordinate;//The integer x y z coordination
        Eigen::Vector3d w_;//Weight that is used in particle to grid function
        Cell() {
            //We should also consider the cell itself for coding simplicity
            this-> v_.setZero();
            this-> tmpV_.setZero();
            this-> w_.setZero();
            this-> pressure = 0;
            this-> type=AIR;
        }
    };

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
                          int num_container_cell_x, int num_container_cell_y){
            this-> container_x_max = container_x_max;
            this-> container_x_min = container_x_min;
            this-> container_y_max = container_y_max;
            this-> container_y_min = container_y_min;

            this-> Nx = num_container_cell_x + 2;//Solid cells around the boundary
            this-> Ny = num_container_cell_y + 2;//Solid cells around the boundary
            this-> N = Nx * Ny;

            this->cell_size_x = (container_x_max - container_x_min) / num_container_cell_x;
            this->cell_size_y = (container_y_max - container_y_min) / num_container_cell_y;

            this-> grid_x_max = container_x_max + cell_size_x;
            this-> grid_x_min = container_x_min - cell_size_x;
            this-> grid_y_max = container_y_max + cell_size_y;
            this-> grid_y_min = container_y_min - cell_size_y;

            this->cells.resize(this->N);

            for(int i = 0; i < cells.size(); i++){
                cells[i] = Cell();
                cells[i].coordinate = this->id_to_coordinate(i);
                if(cells[i].coordinate.x() == 0 || cells[i].coordinate.x() == Nx - 1){
                   cells[i].type = SOLID;
                }
                if(cells[i].coordinate.y() == 0 || cells[i].coordinate.y() == Ny - 1){
                    cells[i].type = SOLID;
                }
                //Add the neighbors
                //Add 4 neighbors in 2d case and 6 neighbors for 3d case
                //x direction neighbor
                for(int nbr_x_idx = std::max(0, cells[i].coordinate.x() - 1);
                    nbr_x_idx <= std::min(Nx - 1, cells[i].coordinate.x() + 1); nbr_x_idx++) {
                    if(nbr_x_idx == cells[i].coordinate.x()){//Skipping the cell itself
                        continue;
                    }
                    cells[i].neighbours.push_back(Eigen::Vector3i(nbr_x_idx,
                                                                  cells[i].coordinate.y(),
                                                                  cells[i].coordinate.z()));
                }
                //y direction neighbor
                for(int nbr_y_idx = std::max(0, cells[i].coordinate.y() - 1);
                    nbr_y_idx <= std::min(Ny - 1, cells[i].coordinate.y() + 1); nbr_y_idx++) {
                    if(nbr_y_idx == cells[i].coordinate.y()){//Skipping the cell itself
                        continue;
                    }
                    //In 3d the 3rd one should also be cells[i].coordinate.z() instead of zero
                    cells[i].neighbours.push_back(Eigen::Vector3i(cells[i].coordinate.x(),
                                                                  nbr_y_idx,
                                                                  cells[i].coordinate.z()));
                }
                //z direction neighbor

                //Adding the cell itself to the neighbors
                cells[i].neighbours.push_back(cells[i].coordinate);

            }
            cell_labels_viewer_idx = Visualize::viewer().append_mesh();
        }

        /*
         * Show the grid inside the viewer window
         */
        void show_grid(){
            // Horizontal grid lines
            for (unsigned y=0; y <= Ny + 1; y++){
                Visualize::viewer().data().add_edges(
                        Eigen::RowVector3d(grid_x_min,
                                           (y - 1) * cell_size_y + grid_y_min,
                                           0),
                        Eigen::RowVector3d(grid_x_max,
                                           (y - 1) * cell_size_y + grid_y_min,
                                           0),
                        Eigen::RowVector3d(0,255,0)
                );
            }
            // Vertical grid lines
            for (unsigned x=0; x <= Nx + 1; x++){
                Visualize::viewer().data().add_edges(
                        Eigen::RowVector3d((x - 1) * cell_size_x + grid_x_min,
                                           grid_y_min,
                                           0),
                        Eigen::RowVector3d((x - 1) * cell_size_x + grid_x_min,
                                           grid_y_max,
                                           0),
                        Eigen::RowVector3d(0,255,0)
                );
            }
        }
        // Functions to work with MAC variables
        /*
         * @breif Convert id of a cell into a 2d id
         * for example id 3 in a 2 by 2 matrix is (1,1) (zero based indexing)
         * Also, it is row vise
         */
        inline Eigen::Vector3i id_to_coordinate(int id){
            Eigen::Vector3i id2D;
            id2D(0) = floor(id / this->Ny);
            id2D(1) = floor(id % this->Ny);
            id2D(2)=0;
            return id2D;
        }

        /*
         * @brief Convert id of a cell into a 2d id
         * for example id 3 in a 2 by 2 matrix is (1,1) (zero based indexing)
         * Also, it is row vise
         */
        inline int coordinate_to_id(Eigen::Vector3i id){
            return id(0) * Ny + id(1);
        }


        /*
         * This function will return the cells real coordination
         * inside the grid
         */
        inline Eigen::Vector3d get_cell_real_coordinate(Eigen::Vector3i coordinate){
            Eigen::Vector3d real_coordination;
            real_coordination(0) = coordinate.x() * cell_size_x + grid_x_min;
            real_coordination(1) = coordinate.y() * cell_size_y + grid_y_min;
            real_coordination(2) == 0;
            assert(real_coordination.x() < grid_x_max);
            assert(real_coordination.y() < grid_y_max);
            assert(real_coordination.x() >= grid_x_min);
            assert(real_coordination.y() >= grid_y_min);
            return real_coordination;

        }
        /*
         * @brief return a pointer to the desire cell.
         * For the sparse implementation, this function should change
         * @param Cell coordination (integer coordination)
         */
        inline Cell& get_cell(Eigen::Vector3i coordination) {
            return cells[coordinate_to_id(coordination)];
        }

        /*
         * @brief given a particle, this function return the cell
         * coordination
         */
        inline Eigen::Vector3i get_particle_cell_coordinate(Eigen::Vector3d prtcl_coordinate){
            Eigen::Vector3i cell_coordinate;
            cell_coordinate(0) = floor((prtcl_coordinate(0) - grid_x_min) / cell_size_x);
            cell_coordinate(1) = floor((prtcl_coordinate(1) - grid_y_min) / cell_size_y);
            cell_coordinate(2) = 0;

            assert(cell_coordinate.x() < Nx);
            assert(cell_coordinate.x() >= 0);
            assert(cell_coordinate.y() < Ny);
            assert(cell_coordinate.y() >= 0);

            return cell_coordinate;
        }
        /*
         * @brief This function shows the cell type on the grid
         */
        void show_cell_type(){
            //Adding type labels
            Visualize::viewer().data_list[cell_labels_viewer_idx].clear_labels();
            for(int i = 0; i < cells.size(); i++){
                std::stringstream l1;
                if(cells[i].type == SOLID)
                    l1 << "S";
                if(cells[i].type == AIR)
                    l1 << "A";
                if(cells[i].type == FLUID)
                    l1 << "F";
                Eigen::Vector3d offset(0.5 * cell_size_x, 0.5 * cell_size_y, 0);
                Visualize::viewer().data_list[cell_labels_viewer_idx].add_label(
                        get_cell_real_coordinate(cells[i].coordinate) + offset, l1.str());
            }
        }

        /*
         * Variables wrapper function
         */
        /*
         * @brief number of grid cells in x direction
         */
        inline int get_Nx() {return Nx;}
        /*
         * @brief number of grid cells in y direction
         */
        inline int get_Ny() {return Ny;}
        /*
         * @brief number of grid cells
         */
        inline int get_N() {return N;}
        /*
         * @brief return the cell size in x direction
         */
        inline double grid_x_size() {return grid_x_max - grid_x_min;}
        /*
         * @brief return the cell size in y direction
         */
        inline double grid_y_size() {return grid_y_max - grid_y_min;}
        /*
         * @brief container boundary size wrapper
         */
        inline double bound_y_max() {return container_y_min;}
        inline double bound_y_min() {return container_y_max;}
        inline double bound_x_max() {return container_x_min;}
        inline double bound_x_min() {return container_x_max;}
        /*
         * @brief get a reference from the cell
         */
        inline std::vector<Cell>& get_cells(){return cells;}
        /*
         * @brief check whether a cell reside
         * in the boundary (in case a particle goes outside of the region
         */
        inline bool is_valid_cell(Eigen::Vector3i& coord){
            if(coord.x() < Nx && coord.x() >= 0 && coord.y() < Ny && coord.y() >= 0)
                return true;
        }
        /*
         * @brief get cell size (for now it returns the size in x direction
         * assuming that the size in y direction is the same
         *
         */
        inline double get_cell_size() {return cell_size_x;}

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

        double grid_x_min;
        double grid_x_max;
        double grid_y_min;
        double grid_y_max;

        int cell_labels_viewer_idx;

        std::vector<Cell> cells;

    };
}
#endif //INC_2D_GRID_H
