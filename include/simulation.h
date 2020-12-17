//
// Created by behrooz on 12/15/20.
//

#ifndef SIMULATION_H
#define SIMULATION_H

#include <visualization.h>
#include <grid.h>
//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>


namespace simulation{
    template <class EigenVector>
    class Sim {
    public:
        Sim() {
          this->macGrid(3,Eigen::Vector2d::Identity());
        }
        ~Sim() = default;
        virtual void simulate_one_step(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, Eigen::VectorXd& forces,
                                       grid::GridDense2D& grid) {};
    private:
      const grid::MacGrid<EigenVector> macGrid;
    };

    template <class EigenVector>
    class GeneralSimulation : public Sim<EigenVector> {
    public:
        /*
         * @brief we use this constructor to work with global variables
         */
        GeneralSimulation() = default;

        ~GeneralSimulation() = default;
        /*
         * @brief This is for initializing the fluid based variables plus gravity
         * @param dt: the initial time step
         * @param g: the gravity acceleration
         * @param density: the density of the fluid
         * @param geometry: The vector of paired (V,F) for obstacles in the scene
         * @param geometry_id: The pair vector of (id, description) of the id
         * of the obstacle in the scene with its definition
         */
        void init_simulation(double dt, EigenVector g, double density, double viscosity,
                             std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> geometry,
                             std::vector<std::pair<int, std::string>> geometry_id){
            this->dt = dt;
            this->gravity = g;
            this->density = density;
            this->vis = viscosity;
        }
        /*
         * @brief This function simulates the fluid from time t^{n} to time t^{n+1}
         * @param q: the general coordination (N*3)
         * @param qdot: the general velocity (N*3)
         */
        void simulate_one_step(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, Eigen::VectorXd& forces,
                               grid::GridDense2D& grid) override {
          this->advect(q, qdot);
          for (int prtc = 0; prtc < q.rows(); prtc++){
            Eigen::Vector3d prtc_q = q.row(prtc);
            Eigen::Vector3d prtc_qdot = qdot.row(prtc);

            // 1) Advect
            prtc_q += prtc_qdot * dt;

            // 2) Update grid
            this->updateGrid();

            // 3a) Adject velocity field
            this->adjectVelocityField();

            // Apply external forces to cells that border fluid cells
            prtc_qdot += gravity * dt;
            q.row(prtc) = prtc_q;
            qdot.row(prtc) = prtc_qdot;

            // Apply viscosity term

            // Calculate pressure

            // Apply pressure

            // Extrapolate fluid velocities into surrounding cells
            for (auto& cell : this->macGrid) {
              cell.layer = cell.type == grid::CellTypes::FLUID ? 0 : -1;
            }
            for (int i = 1; i < std::max(2, kfcl); i++) {
              for (auto& cell : this->macGrid) {
                if (cell.layer != -1)
                  continue;
                // if cell has neighbour with layer == i-1
              }
            }

            // Set the velocities that point to solid cells
            for (auto& cell : this->macGrid) {
            }

            // Move the marker particles through the fluid
            for (auto& p : this->particles) {
              // Determine if frame should be rendered
              if (this->time % this->fps) {
                p.p = p.p + p.v * this->dt;
              } else {
                // RK4 Method
              }
            }
          }
        }

      /*
       * Perform global update to macgrid prior to calculations
       */
      void updateGrid() {
        // Grid update
        auto particles = this->macGrid.getParticles();
        auto numElm = particles.size() / this->macGrid.numParticles;
        for (int i = 0; i < this->macGrid.numParticles; i++) {
          // Get cell in which particle exists
          Eigen::Ref<EigenVector> p = particles.segment<numElm>(i*numElm);
          auto* cell = this->macGrid.getCellFromParticle(p);

          // Check if cell is not within defined bounds
          if (cell == nullptr) {
            if (this->cellIsWithinSimBounds()) {
              typename grid::MacGrid<EigenVector>::Cell newCell = { grid::CellTypes::FLUID, 0 };
              this->macGrid[p] = cell;
            }
          } else if (this->cellIsInSolidBounds(cell)) {
            cell->type = grid::CellTypes::FLUID;
            cell->layer = 0;
          }
        }

        // Create a buffer zone around the fluid
        for (int i = 1; i < std::max(2, std::ceil(kcfl)); i++) {
          for (auto& cell : this->macGrid) {
            if (cell.layer != i-1) {
              continue;
            }
            // Loop through neighbours of cell
            auto neighbours = this->macGrid.generateNeighboursFromCell(cell);

            // TODO: Update pointer connections
            for (auto& n : neighbours) {
              auto* nc = this->macGrid.getCellFromParticle(n);
              if (nc != nullptr) {
                if (nc->layer == -1 && nc->type != grid::CellTypes::SOLID) {
                  nc->type = grid::CellTypes::AIR;
                  nc->layer = i;
                }
              } else {
                this->macGrid[nc] = grid::MacGrid<EigenVector>::Cell{
                    i,
                    this->cellIsWithinSimBounds() ? grid::CellTypes::AIR : grid::CellTypes::SOLID
                };
              }
            }
          }
        }
        // TODO: Delete any cell with layer == -1
      }

      void adjectVelocityField() {
        for (auto& cell : this->macGrid) {

        }
      }

      EigenVector traceParticle(EigenVector& e) {
        auto v = this->macGrid.interpolateVelocity(e);
        auto v2 = this->macGrid.interpolateVelocity(e + v*0.5);

        return e + v2 * this->dt;
      }

      bool cellIsInSolidBounds(){}
      bool cellIsWithinSimBounds(){}
    private:
      EigenVector gravity;
      double dt;
      double density;
      double vis;

      // Simulation helper functions
      inline void advect(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot){
      }
    };

//Simulation objects
GeneralSimulation<Eigen::Vector2d> water_sim;
}

#endif //A6_RIGID_BODIES_CONTACT_SIMULATION_H
