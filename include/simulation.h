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


    class GeneralSimulation{
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
                             std::vector<std::pair<int, std::string>> geometry_id,
                             Eigen::MatrixXd& q, Eigen::MatrixXd& qdot){
            this->dt = dt;
            this->gravity = g;
            this->density = density;
            this->vis = viscosity;
            //Get the rigid coordination
            for(int rig_id = 0; rig_id < geometry_id.size(); rig_id++){
                if(geometry_id[rig_id].second == "boundary"){
                    boundary_V = geometry[rig_id].first;
                    boundary_F = geometry[rig_id].second;

                }
            }
            //getting the boundaty sizes
            Eigen::Vector3d m = boundary_V.colwise().minCoeff();
            Eigen::Vector3d M = boundary_V.colwise().maxCoeff();
            x_size_max = M(0);
            y_size_max = M(1);
            x_size_min = m(0);
            y_size_min = m(1);
            //Init the Grid
            //Mark the particles cell
            int num_container_cell_x = 10;
            int num_container_cell_y = 10;
            grid2d.init_2d_grid(x_size_min, x_size_max,
                                y_size_min, y_size_max,
                                num_container_cell_x, num_container_cell_y);

            // Create The Fluid inside the box (a rectangle shape points)
            num_particles = 100;
            // Start to create the particles in the top left of the square
            double fluid_x_range = (x_size_max - x_size_min) / 8.0;
            double fluid_y_range = (y_size_max - y_size_min) / 2.0;
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
            particles_viewer_id = Visualize::viewer().append_mesh();
            Visualize::viewer().data_list[particles_viewer_id].set_points(q, Eigen::RowVector3d(0,0,128));
            Visualize::viewer().data().point_size = 5;
            // Updating grid based on the marker
//
//            for(int x = 1; x < num_cell_x / 2; x++){
//                for(int y = 1; y < num_cell_y / 2; y++) {
//
//                }
//            }

        }
        /*
         * @brief This function simulates the fluid from time t^{n} to time t^{n+1}
         * So we are doing the following steps based on the Fluid flow for the rest of us: Tutorial of the marker and cell method in computer graphics
         * and using the book tutorial and the ETHz course notes and the lecture video
         * 1. compute the velocity field (particles to grid - classifying the cells - saving the field for FLIP step)
         * 2. apply external forces
         * 3. apply viscosity
         * 4. Enforce boundary conditions
         * 5. compute and apply pressure gradients
         * 6. update particles velocities
         * 7. get the dt correct for sparse grid
         * 8. move the particles and check for boundaries
         * @param q: the general coordination (N*3)
         * @param qdot: the general velocity (N*3)
         * @param forces: is the external forces
         */
        void simulate_one_step(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, Eigen::VectorXd& forces) {
            // Map the particles velocity into the cells
            this->particles_to_cell(q, qdot);

            // Saving the V for FLIP step and apply Forces
            for(auto &iter: grid2d.get_cells()){
                iter.tmpV_ = iter.v_;
                iter.v_ += gravity * dt;
            }

            // Apply viscosity
//            this->apply_viscosity();

            // Apply pressure
            this->apply_pressure();

            // Enforce boundary conditions
            this->enforce_boundaries();

            // Update particles
            this->cells_to_particles(q, qdot);

            // Calculate the next dt based on the fluids note paper
            // Doesn't matter for now

            // Move the particles Advect
            for (int prtc = 0; prtc < q.rows(); prtc++){
                Eigen::Vector3d prtc_q = q.row(prtc);
                Eigen::Vector3d prtc_qdot = qdot.row(prtc);
                prtc_q += prtc_qdot * dt;
                q.row(prtc) = prtc_q;
            }
        }

        /*
         * Enforce boundaries based on where solid cells are
         */
        void enforce_boundaries() {
            auto cells = this->grid2d.get_cells();
            for (auto& c : cells) {
                for (auto& nCoord : c.neighbours) {
                    auto nc = this->grid2d.get_cell(nCoord);

                    // Determine if neighbouring cell is solid
                    if (nc.type == grid::SOLID) {

                        // Determine direction to cell
                        auto dir = (nc.coordinate - c.coordinate).cast<double>();

                        /*
                         * If sum of the dot-product is
                         * greater than zero, the velocity vector
                         * points into the solid cell.
                         */
                        double sum = c.v_.dot(dir);
                        if (sum > 0) {
                            c.v_ -= (dir * sum);
                        }
                    }
                }
            }
        }

        /*
         * Update cells velocity based on viscosity
         */
        void apply_viscosity() {
            // Obtain cells
            auto cells = this->grid2d.get_cells();

            // Setup update vector
            Eigen::VectorXd vNew;
            vNew.resize(cells.size() * 3);

            for (int i = 0; i < cells.size(); i++) {
                auto cell = cells[i];
                vNew.segment<3>(i*3) = this->dt * this->vis * this->calculateLaplacian(cell);
            }
        }

        /*
         * Calculate laplacian for given cell
         */
        Eigen::Vector3d calculateLaplacian(grid::Cell& c) {
            Eigen::Vector3d laplacian;

            // Obtain all non solid neighbours
            for (int dim = 0; dim < 3; dim++) {
                // Determine if dim borders fluid cells

                // TODO: Determine if 1/-1 are adjacent cells
                auto key = c.coordinate;
                bool isBorderingFluidCell = false;
                for (int k = 0; k < 2; k++) {
                    key(dim) += k == 0 ? (-1) * 1 : 1;
                    auto nc = this->grid2d.get_cell(key);
                    if (nc.type == grid::CellTypes::FLUID) {
                        isBorderingFluidCell = true;
                        break;
                    }
                    key(dim) += k == 1 ? (-1) * 1 : 1;
                }

                // Calculate the laplacian
                if (isBorderingFluidCell) {
                    for (auto& nCoords : c.neighbours) {
                        auto nc = this->grid2d.get_cell(nCoords);

                        if (nc.type == grid::CellTypes::FLUID) {
                            laplacian(dim) += nc.v_(dim) - c.v_(dim);
                        }
                    }
                }
            }
            return laplacian;
        }

        /*
         * Update cells velocity based on pressure
         */
        void apply_pressure() {
            // Setup sparse matrix triplets
            std::vector<Eigen::Triplet<double>> triplets;

            // Get cells we are working with
            auto cells = this->grid2d.get_cells();

            // Setup rhs vector
            Eigen::VectorXd b;
            b.resize(cells.size());

            // We use the unique cell identifiers to build matrix
            for (int i = 0; i < cells.size(); i++) {
                auto cell = cells[i];

                int nSolid = 0;
                double aCells = 0.;
                for (auto& n : cell.neighbours) {
                    auto nc = grid2d.get_cell(n);

                    if (nc.type != grid::CellTypes::SOLID) {
                        nSolid++;
                    }
                    if (nc.type == grid::CellTypes::FLUID) {
                        triplets.push_back(Eigen::Triplet<double>{
                            this->grid2d.coordinate_to_id(cell.coordinate),
                            this->grid2d.coordinate_to_id(nc.coordinate),
                            1.
                        });
                    }
                    if (nc.type == grid::CellTypes::AIR) {
                        aCells += 1.;
                    }
                }
                // Assemble A
                triplets.push_back(Eigen::Triplet<double>{
                    this->grid2d.coordinate_to_id(cell.coordinate),
                    this->grid2d.coordinate_to_id(cell.coordinate),
                    (-1.) * nSolid
                });

                // Calculate entry for vector b
                // TODO: replace with actual height
                b(i) =  (this->density * /**/ 1) *
                        (this->calculate_divergence(cell) / this->dt) -
                        (aCells * this->atm);
            }

            // Create solution matrix A
            Eigen::SparseMatrix<double> A;
            A.resize(this->grid2d.get_N(), this->grid2d.get_N());
            A.setFromTriplets(triplets.begin(), triplets.end());

            // Solve system
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;

            cg.compute(A);

            if (cg.info() != 0) {
                std::runtime_error("Matrix solver did not work...");
            }

            Eigen::VectorXd x = cg.solve(b);

            // Apply pressure update
            for (int i = 0; i < cells.size(); i++) {
                auto cell = cells[i];
                Eigen::Vector3d border = {1.,1.,1.};
                for (auto nCoord : cell.neighbours) {
                    auto nc = this->grid2d.get_cell(nCoord);
                    if (nc.type == grid::CellTypes::SOLID) {
                        auto diff = cell.coordinate - nCoord;
                        for (int k = 0; k < diff.size(); k++) {
                            if (diff(k) != 0) {
                                border(k) = 0;
                            }
                        }
                    }
                }
                cells[i].v_ = cells[i].v_ - border * this->dt / (this->density * /* height */ 1) * x(i);
            }
        }

        /*
         * Calculates divergence at given cell
         */
        double calculate_divergence(grid::Cell& c) {
            auto key = c.coordinate;

            /*
             * TODO: Update using neighbour code
             */
            double divergence = 0.;
            for (int i = 0; i < key.size(); i++) {
                // Look one cell ahead
                key(i) += 1;

                if (key(i) < this->grid2d.get_Nx()) {
                    auto nc = this->grid2d.get_cell(key);
                    if (nc.type != grid::CellTypes::SOLID) {
                        divergence += nc.v_(i);
                    }
                }

                // Look one cell behind
                key(i) -= 2;

                if (key(i) >= 0) {
                    auto nc = this->grid2d.get_cell(key);
                    if (nc.type != grid::CellTypes::SOLID) {
                        divergence -= c.v_(i);
                    }
                }

                // Reset key
                key(i) += 1;
            }
            return divergence;
        }

        /*
         * Update cells based on particles position
         */
        void particles_to_cell(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot){
            //Preparing the cells for the simulation in this step
            for(auto &iter: this->grid2d.get_cells()){
                iter.v_.setZero();
                iter.w_.setZero();
                iter.pressure = 0;
                iter.type = grid::AIR;
                if(iter.coordinate.x() == 0 || iter.coordinate.x() == grid2d.get_Nx() - 1){
                    iter.type = grid::SOLID;
                }
                if(iter.coordinate.y() == 0 || iter.coordinate.y() == grid2d.get_Ny() - 1){
                    iter.type = grid::SOLID;
                }
            }
            //For each particle
            for(int prtcl = 0; prtcl < q.rows(); prtcl++){
                auto cell = grid2d.get_cell(grid2d.get_particle_cell_coordinate(q.row(prtcl)));
                if(cell.type == grid::SOLID) {std::cout << "The particle is not in the right location" << std::endl;}
                Eigen::Vector3d prtcl_q = q.row(prtcl);
                Eigen::Vector3d prtcl_qdot = qdot.row(prtcl);
                //apply velocity to all the neighboring cells (radius 1 for now)
                for(auto &iter: cell.neighbours){
                    //We still have to consider the solid cells that have velocity components
                    //Toward the boundary.
                    apply_velocity_to_cell(grid2d.get_cell(iter), prtcl_q, prtcl_qdot);
                }
            }

            // Normalize the velocity
            // can have #pragma
            for(int id = 0; id < grid2d.get_N(); id++){
                auto cell = grid2d.get_cell(grid2d.id_to_coordinate(id));
                if(cell.type == grid::FLUID){
                    if(cell.w_.x() != 0) cell.v_.x() = cell.v_.x() / cell.w_.x();
                    if(cell.w_.y() != 0) cell.v_.y() = cell.v_.y() / cell.w_.y();
//                    cell.v_.z() = cell.v_.z() / cell.w_.z();
                }
            }

            //This also can have pragma
            for(int id = 0; id < grid2d.get_N(); id++){
                auto cell = grid2d.get_cell(grid2d.id_to_coordinate(id));
                int cnt = 0;
                if(cell.type == grid::AIR){
                    for(auto& iter: cell.neighbours){
                        auto nbr = grid2d.get_cell(iter);
                        if(nbr.type == grid::FLUID){//Extrapolate
                            cell.v_ = nbr.v_;
                            cnt++;
                        }
                    }
                }
                cell.v_ = cell.v_ / cnt;
            }
        }

        /*
         * @brief given a particle and a cell, this function applies the velocity
         * to the cell. Note that for 3d this function should change a bit
         * @param cell: The cell that we want to apply force to
         * @param prtcl_q is the position of the particle
         * @param prtcle_qdot is the velocity of the particle
         */
        inline void apply_velocity_to_cell(grid::Cell& cell, Eigen::Vector3d& prtcl_q, Eigen::Vector3d& prtcl_qdot){
            //apply the velocity to Vx, Vy, Vz separately
            //Applying to Vx
            //Now we only extrapolate on fluid cells
            cell.type = grid::FLUID;
            Eigen::Vector3d u;
            //The cell coordinate is the bottom left vertex of the cell
            u.x() = (cell.coordinate.x()) * grid2d.grid_x_size();
            u.y() = (cell.coordinate.y() + 0.5) * grid2d.grid_y_size();
            u.z() = 0;
            double weight = sph_weight(prtcl_q, u, 2 * grid2d.grid_x_size());
            if(weight != 0){
                cell.v_.x() += weight * prtcl_qdot.x();
                cell.w_.x() += weight;
            }
            //Applying to Vy
            u.x() = (cell.coordinate.x() + 0.5) * grid2d.grid_x_size();
            u.y() = (cell.coordinate.y()) * grid2d.grid_y_size();
            u.z() = 0;
            weight = sph_weight(prtcl_q, u, 2 * grid2d.grid_y_size());
            if(weight != 0){
                cell.v_.y() += weight * prtcl_qdot.y();
                cell.w_.y() += weight;
            }
            //Applying to Vz
        }


        /*
         * @brief it is poly6 weight function according to
         * the ETHz slides (one of the SPH kernels)
         */
        inline double sph_weight(Eigen::Vector3d& prtcl_pos, Eigen::Vector3d& v_pos, double h){
            double dist = (prtcl_pos - v_pos).norm();
            if(dist < h){
                double diff = std::pow((std::pow(h, 2) - std::pow(dist, 2)), 3);
                return (315/(64 * M_PI * std::pow(h, 9))) * diff;
            } else {
                return 0;
            }
        }

        inline void cells_to_particles(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot){
            // For each particle
            for (int prtcl = 0; prtcl < q.rows(); prtcl++){
                auto cell = grid2d.get_cell(grid2d.get_particle_cell_coordinate(q.row(prtcl)));
                if (cell.type == grid::SOLID) {std::cout << "The particle is not in the right location" << std::endl;}
                Eigen::Vector3d prtcl_q = q.row(prtcl);
                // w(0) is the weight of Vx1. w(0) use the distance from Vx2
                // w(1) is the weight of Vx2. w(1) use the distance from Vx1

                // w(2) is the weight of Vy1. w(2) use the distance from Vy2
                // w(3) is the weight of Vy2. w(3) use the distance from Vy1
                Eigen::Vector4d w;
                Eigen::Vector4d v;
                Eigen::Vector4d tmpv;
                v.setZero();
                tmpv.setZero();
                w.setZero();

                w(1) = prtcl_q.x() - cell.coordinate.x() * grid2d.get_cell_size();
                w(3) = prtcl_q.y() - cell.coordinate.y() * grid2d.get_cell_size();

                v(0) = cell.v_(0);
                tmpv(0) = cell.tmpV_(0);

                v(2) = cell.v_(1);
                tmpv(2) = cell.tmpV_(1);

                // Apply velocity to particle
                for (auto &iter: cell.neighbours){
                    // bilinear interpolation between the particle cell and the cell in the right
                    if (grid2d.get_cell(iter).coordinate.x() == cell.coordinate.x() + 1){
                        w(0) = prtcl_q.x() - grid2d.get_cell(iter).coordinate.x() * grid2d.get_cell_size();
                        v(1) = grid2d.get_cell(iter).v_(0);
                        tmpv(1) = grid2d.get_cell(iter).tmpV_(0);
                    }
                    if (grid2d.get_cell(iter).coordinate.y() == cell.coordinate.y() + 1){
                        w(2) = prtcl_q.y() - grid2d.get_cell(iter).coordinate.y() * grid2d.get_cell_size();
                        v(3) = grid2d.get_cell(iter).v_(1);
                        tmpv(3) = grid2d.get_cell(iter).tmpV_(1);
                    }
                }

                // FLIP operation
                qdot.row(prtcl).x() += w(0) * (v(0) - tmpv(0)) + w(1) * (v(1) - tmpv(1));
                qdot.row(prtcl).y() += w(2) * (v(2) - tmpv(2)) + w(3) * (v(3) - tmpv(3));
                qdot.row(prtcl).z() = 0;
            }
        }

        //============ Wrapper functions ==============
        /*
         * Get the layer id for the particles marker in the viewer
         */
        int get_prtcl_veiwer_id() {return particles_viewer_id;}
        /*
         * Get the grid object
         */
        inline grid::GridDense2D& get_grid(){ return grid2d; }

    private:
        grid::GridDense2D grid2d;
        Eigen::Vector3d gravity;
        double dt;
        double density;
        double vis;
        Eigen::MatrixXd boundary_V;
        Eigen::MatrixXi boundary_F;
        double x_size_min;
        double x_size_max;
        double y_size_min;
        double y_size_max;
        int num_particles;
        int particles_viewer_id;
        double atm;
    };
}

#endif //A6_RIGID_BODIES_CONTACT_SIMULATION_H
