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
    class Sim{
    public:
        Sim() = default;
        ~Sim() = default;
        virtual void simulate_one_step(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, Eigen::VectorXd& forces,
                                       grid::GridDense2D& grid) {};
    private:
    };

    class GeneralSimulation: public Sim{
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
        void init_simulation(double dt, Eigen::Vector3d g, double density, double viscosity,
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
            for(int prtc = 0; prtc < q.rows(); prtc++){
                Eigen::Vector3d prtc_q = q.row(prtc);
                Eigen::Vector3d prtc_qdot = qdot.row(prtc);
                // Advect
                prtc_q += prtc_qdot * dt;
                // External Forces
                prtc_qdot += gravity * dt;
                q.row(prtc) = prtc_q;
                qdot.row(prtc) = prtc_qdot;
            }
        }
    private:
        Eigen::Vector3d gravity;
        double dt;
        double density;
        double vis;
        // Simulation helper functions
        inline void advect(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot){

        }
    };

    //Simulation objects
    GeneralSimulation water_sim;

}

#endif //A6_RIGID_BODIES_CONTACT_SIMULATION_H
