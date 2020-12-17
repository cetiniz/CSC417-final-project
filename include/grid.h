//
// Created by behrooz on 12/16/20.
//


#ifndef INC_2D_GRID_H
#define INC_2D_GRID_H

#include <fluid_utils.h>
#include <visualization.h>

namespace grid{
    /*
     * Different types for the cells inside a macgrid
     */
    enum CellTypes {
      SOLID,
      FLUID,
      AIR
    };


    /*
     * Data structure to represent macgrid of fluid simulation.
     *
     * EigenVector -> [Eigen::Vector2d || Eigen::Vector3d]
     */
    template <class EigenVector>
    class MacGrid {
      /*
       * Struct containing data for a given cell of the macgrid
       */
      struct Cell {
        CellTypes type;
        int layer;
        std::vector<Cell*> neighbours;
        EigenVector v;
        EigenVector tmpV;

        Cell() {
          this->neighbours.reserve(8);
        }
      };

      static const bool eq(const EigenVector& lhs, const EigenVector& rhs){ return lhs == rhs; };
      static const std::size_t hash (const EigenVector & n){ return (std::size_t)(n*Eigen::Vector3d{541 + 79 + 31});};
      typedef std::unordered_map<EigenVector, Cell*, decltype(hash), decltype(eq)> SparseGrid;

    public:
      MacGrid(int h, EigenVector o) :
          height(h),
          origin(o)
      {}

      typename SparseGrid::iterator begin() {
        return this->sparseGrid.begin();
      }

      typename SparseGrid::iterator end() {
        return this->sparseGrid.end();
      }

      bool containsParticle(EigenVector& e) {
        return this->sparseGrid.count(e) > 0;
      }

      bool containsCell(Cell* c) {
        return this->sparseGrid.count(c.origin);
      }

      Cell* getCellFromParticle(EigenVector& e) {

      }

      EigenVector interpolateVelocity(EigenVector& e) {
        EigenVector ne;
        for (int i = 0; i < e.size(); i++) {
          ne = (e / this->height) - (EigenVector::Identity() * 0.5);
          ne(i) += 0.5;
          ne(i) = this->_calculateInterpretedVelocity(ne,i);
        }
      }

      double _calculateInterpretedVelocity(EigenVector ne, int index) {
        if (ne.size() == 3) {
          int i = std::floor(ne(0));
          int j = std::floor(ne(1));
          int k = std::floor(ne(2));

          // @TODO - finish interpreted velocity
          return 0.0;
        } else if (ne.size() == 2) {

        }
      }

      std::vector<Cell*> generateNeighboursFromCell(Cell* c) {

      }

      Eigen::Ref<Eigen::VectorXd> getParticles() {
        return this->markers;
      }

      Cell* operator[](EigenVector& e) {
        return sparseGrid[e];
      }

    private:
      Eigen::VectorXd markers;
      EigenVector origin;
      double height{};
      SparseGrid sparseGrid;
    };

class Grid {
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
