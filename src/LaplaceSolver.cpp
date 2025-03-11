/* ------------------------------------------------------------------------
*
 * SPDX-License-Identifier: GPL-3-or-later
 * Copyright (C) 2023 - 2024 by Matteo Malvestiti
 *
 * This file is part of ion_propulsion.
 *
 * ion_propulsion is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ion_propulsion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ion_propulsion; see the file COPYING.  If not, see
 * <https://www.gnu.org/licenses/>.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Matteo Malvestiti, Politecnico di Milano, 2024
 *
 */

#include "LaplaceSolver.h"
#include <Constants.h>

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{

    // ------------------------------------------------------
    // Base
    // ------------------------------------------------------
    template <int dim>
    Base<dim>::Base(parallel::distributed::Triangulation<dim> &coarse_grid)
      : triangulation(&coarse_grid)
      , convergence_table(std::make_shared<ConvergenceTable>())
      , refinement_cycle(numbers::invalid_unsigned_int)

      , mpi_communicator(MPI_COMM_WORLD)
      , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
      , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
      , pcout(std::cout, (this_mpi_process == 0))

      , computing_timer(this->mpi_communicator,
                      this->pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
    {}

    template <int dim>
    void Base<dim>::set_refinement_cycle(const unsigned int cycle)
    {
      refinement_cycle = cycle;
    }

    template <int dim>
    void Base<dim>::checkpoint()
    {
      std::cout << "--- Writing checkpoint... ---" << std::endl;

//#if (DEAL_II_VERSION_MAJOR > 9) || (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR >= 6)
      triangulation->save(OUTPUT_PATH+"/"+"checkpoint-mesh"+std::to_string(this->refinement_cycle));
//#else
      //std::ofstream                 checkpoint_file(OUTPUT_PATH+"/"+"checkpoint-mesh"+std::to_string(this->refinement_cycle));
      //boost::archive::text_oarchive archive(checkpoint_file);
      //triangulation->save(archive,1);
//#endif
    }

    template <int dim>
    void Solver<dim>::restart()
      {
        {
          std::ifstream checkpoint_file(OUTPUT_PATH+"/"+"checkpoint_ion_propulsion");
          AssertThrow(checkpoint_file,
                      ExcMessage(
                        "Could not read from the <checkpoint_ion_propulsion> file."));

          boost::archive::text_iarchive archive(checkpoint_file);
          archive >> *this;
        }

      this->triangulation->load(OUTPUT_PATH+"/"+"checkpoint");
      //particle_handler.deserialize();

      dof_handler.reinit(*this->triangulation);

      GridOut grid_out;
      GridOutFlags::Msh msh_flags(true, true);
      grid_out.set_flags(msh_flags);
      grid_out.write_msh(*this->triangulation, OUTPUT_PATH+"/re-imported_mesh.msh");
    }


    // ------------------------------------------------------
    // Solver
    // ------------------------------------------------------

    template <int dim>
    Solver<dim>::Solver(parallel::distributed::Triangulation<dim> &       triangulation,
                        const FiniteElement<dim> & fe,
                        const Quadrature<dim> &    quadrature,
                        const Quadrature<dim - 1> &face_quadrature,
                        const Function<dim> &      boundary_values,
                        const unsigned degree)
      : Base<dim>(triangulation)
      , fe(&fe)
      , quadrature(&quadrature)
      , face_quadrature(&face_quadrature)
      , dof_handler(triangulation)
      , boundary_values(&boundary_values)
      , mapping(degree)
    {}


    template <int dim>
    Solver<dim>::~Solver()
    {
      dof_handler.clear();
    }

    template <int dim>
    void Solver<dim>::update_convergence_table() {
      if (this->this_mpi_process == 0) {
        this->convergence_table->add_value("cycle", this->refinement_cycle);
        //this->convergence_table->add_value("cells", this->triangulation->n_active_cells());
        this->convergence_table->add_value("DoFs", this->dof_handler.n_dofs());

        CSVLogger& logger = CSVLogger::getInstance();
        logger.addColumn("cycle", std::to_string(this->refinement_cycle));
        //logger.addColumn("cells", std::to_string(this->triangulation->n_active_cells()));
        logger.addColumn("DoFs", std::to_string(this->dof_handler.n_dofs()));
      }
    }

    template <int dim>
    void Solver<dim>::solve_problem()
    {
      // SETUP scope
      {
        TimerOutput::Scope t(this->computing_timer, "setup");

        dof_handler.distribute_dofs(*fe);
        this->locally_owned_dofs = dof_handler.locally_owned_dofs();
        this->locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

        locally_relevant_solution.reinit(this->locally_owned_dofs,
                                         this->locally_relevant_dofs,
                                         this->mpi_communicator);
        completely_distributed_solution.reinit(this->locally_owned_dofs,
                                                this->mpi_communicator);


        if(MANUAL_LIFTING_ON)     //TODO
          if (this->refinement_cycle == 0) {
            locally_relevant_Rg_vector.reinit(this->locally_owned_dofs, this->locally_relevant_dofs, this->mpi_communicator);
            construct_Rg_vector();
          }


        linear_system_ptr = std::make_unique<LinearSystem>(
                  dof_handler,
                  this->mpi_communicator,
                  this->locally_owned_dofs,
                  this->locally_relevant_dofs);
      }

      // ASSEMBLE scope
      {
        TimerOutput::Scope t(this->computing_timer, "assembly");
        assemble_linear_system(*linear_system_ptr);
      }

      // SOLVE scope
      {
        TimerOutput::Scope t(this->computing_timer, "solve");
        linear_system_ptr->solve(locally_relevant_solution, this->completely_distributed_solution);

        if(MANUAL_LIFTING_ON)   //TODO
          retrieve_Rg();
      }

    }

    template <int dim>
    void Solver<dim>::postprocess(
      const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      std::pair<std::string, double> possible_pair = postprocessor(dof_handler, locally_relevant_solution,*this->triangulation);
      if (possible_pair.first != "null") {
        Base<dim>::convergence_table->add_value(possible_pair.first, possible_pair.second);
        Base<dim>::convergence_table->set_scientific(possible_pair.first, true);
        CSVLogger::getInstance().addColumn(possible_pair.first, to_string_with_precision(possible_pair.second,15));
      }

    }


    template <int dim>
    unsigned int Solver<dim>::n_dofs() const
    {
      return dof_handler.n_dofs();
    }

    template <int dim>
    void Solver<dim>::assemble_linear_system(LinearSystem &linear_system) {

      linear_system.rhs    = 0;
      linear_system.matrix = 0;
      linear_system.Umatrix = 0;

      // ------------------------------------------------------
      // Matrix assemble
      // ------------------------------------------------------
      FEValues<dim> fe_values(mapping, *fe,
                              *quadrature,
                              update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
      const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
      const unsigned int n_q_points    = quadrature->size();
      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      std::vector<double>  rhs_values;
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned()){
          cell_matrix = 0;
          fe_values.reinit(cell);

          for (unsigned int i = 0; i < dofs_per_cell; ++i){
            for (unsigned int j = 0; j < dofs_per_cell; ++j){
              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point){
                cell_matrix(i, j) += eps_0 * eps_r *
                    fe_values.shape_grad(i, q_point) *
                    fe_values.shape_grad(j, q_point) *
                    fe_values.JxW(q_point);
                }
              }
            }
          cell->get_dof_indices(local_dof_indices);
          linear_system.hanging_node_constraints.distribute_local_to_global(
            cell_matrix,
            local_dof_indices,
            linear_system.matrix);
        }
      linear_system.matrix.compress(VectorOperation::add);

      // ------------------------------------------------------
      // RHS assemble
      // ------------------------------------------------------
      assemble_rhs(linear_system.rhs, linear_system.hanging_node_constraints);
      linear_system.rhs.compress(VectorOperation::add);

      // ------------------------------------------------------
      // BCs
      // ------------------------------------------------------
      std::map<types::global_dof_index, double> boundary_value_map;
      interpolate_boundary_values(boundary_value_map);

      MatrixTools::apply_boundary_values(boundary_value_map,
                                         linear_system.matrix,
                                         completely_distributed_solution,
                                         linear_system.rhs,
                                         false);

    }

    template <int dim>
        Solver<dim>::LinearSystem::LinearSystem(
          const DoFHandler<dim> &dof_handler,
          MPI_Comm & mpi_communicator,
          IndexSet & locally_owned_dofs,
          IndexSet & locally_relevant_dofs)
    {
      hanging_node_constraints.clear();

#if (DEAL_II_VERSION_MAJOR > 9) || (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR >= 6)
      hanging_node_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
#else
      hanging_node_constraints.reinit(locally_owned_dofs);
#endif
      DoFTools::make_hanging_node_constraints(dof_handler, hanging_node_constraints);
      hanging_node_constraints.close();

      DynamicSparsityPattern dsp(locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler,
                                      dsp,
                                      hanging_node_constraints,
                                      /*keep constrained dofs*/ false);
      SparsityTools::distribute_sparsity_pattern(dsp,
                                                 dof_handler.locally_owned_dofs(),
                                                 mpi_communicator,
                                                 locally_relevant_dofs);

      matrix.reinit(locally_owned_dofs,
                    locally_owned_dofs,
                    dsp,
                    mpi_communicator);
      Umatrix.reinit(locally_owned_dofs,
                    locally_owned_dofs,
                    dsp,
                    mpi_communicator);
      rhs.reinit(locally_owned_dofs, mpi_communicator);

    }



    template <int dim>
    void Solver<dim>::LinearSystem::solve(
      PETScWrappers::MPI::Vector &locally_relevant_solution,
      PETScWrappers::MPI::Vector & completely_distributed_solution
      ) const
    {


      SolverControl            solver_control(5000, 1e-12);
      PETScWrappers::SolverCG  cg(solver_control);

      PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
      data.symmetric_operator = true;
      PETScWrappers::PreconditionBoomerAMG preconditioner;
      preconditioner.initialize(matrix, data);

      cg.solve(matrix,
                completely_distributed_solution,
                rhs,
                preconditioner);

      /*PETScWrappers::SparseDirectMUMPS solverMUMPS(solver_control);
      solverMUMPS.solve(matrix, completely_distributed_solution, rhs);*/

      hanging_node_constraints.distribute(completely_distributed_solution);

      locally_relevant_solution = completely_distributed_solution;

      cout<<"Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<std::endl;
    }

    // ------------------------------------------------------
    // PrimalSolver
    // ------------------------------------------------------

    template <int dim>
    PrimalSolver<dim>::PrimalSolver(parallel::distributed::Triangulation<dim> &       triangulation,
                                    const FiniteElement<dim> & fe,
                                    const Quadrature<dim> &    quadrature,
                                    const Quadrature<dim - 1> &face_quadrature,
                                    const Function<dim> &      rhs_function,
                                    const Function<dim> &      boundary_values,
                                    const unsigned degree)
      : Base<dim>(triangulation)
      , Solver<dim>(triangulation,
                    fe,
                    quadrature,
                    face_quadrature,
                    boundary_values,
                    degree)
      , rhs_function(&rhs_function)
    {}

    template <int dim>
    void PrimalSolver<dim>::output_solution()
    {
      TimerOutput::Scope t(this->computing_timer, "output");

      DataOut<dim> data_out;
      data_out.attach_dof_handler(this->dof_handler);
      data_out.add_data_vector(this->locally_relevant_solution, "uh",DataOut<dim, dim>::type_dof_data);


      Vector<float> subdomain(this->triangulation->n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = this->triangulation->locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");

      /*Vector<double> rhs_function_values(this->dof_handler.n_dofs());
      VectorTools::interpolate(this->dof_handler, *this->rhs_function, rhs_function_values);
      data_out.add_data_vector(rhs_function_values, "rhs_function",DataOut<dim, dim>::type_dof_data);

      Vector<double> uex_function_values(this->dof_handler.n_dofs());
      VectorTools::interpolate(this->dof_handler, *this->boundary_values, uex_function_values);
      data_out.add_data_vector(uex_function_values, "u_ex",DataOut<dim, dim>::type_dof_data);

      Vector<double> boundary_ids(this->triangulation->n_active_cells());
      for (const auto &cell : this->triangulation->active_cell_iterators())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            boundary_ids[cell->active_cell_index()] = cell->face(face)->boundary_id();
      data_out.add_data_vector(boundary_ids, "boundary_ids",DataOut<dim, dim>::type_cell_data);*/

      data_out.build_patches(this->mapping, 1,DataOut<dim,dim>::CurvedCellRegion::curved_inner_cells);

      data_out.write_vtu_with_pvtu_record(
        OUTPUT_PATH+"/", "solution-", this->refinement_cycle, this->mpi_communicator, 2);
    }

    template <int dim>
    void PrimalSolver<dim>::assemble_rhs(PETScWrappers::MPI::Vector &rhs, AffineConstraints<double> & constraints) const
    {
      FEValues<dim> fe_values(this->mapping, *this->fe,
                              *this->quadrature,
                              update_values | update_gradients | update_quadrature_points |
                                update_JxW_values);
      const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
      const unsigned int n_q_points    = this->quadrature->size();

      Vector<double>                       cell_rhs(dofs_per_cell);
      std::vector<double>                  rhs_values(n_q_points);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Tensor<1, dim>>          rg_gradients(n_q_points);

      for (const auto &cell : this->dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned()) {
          cell_rhs = 0;

          fe_values.reinit(cell);
          rhs_function->value_list(fe_values.get_quadrature_points(),
                                   rhs_values);
          if(MANUAL_LIFTING_ON)
            fe_values.get_function_gradients(this->locally_relevant_Rg_vector, rg_gradients);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                              rhs_values[q_point] *               // f((x_q)
                              fe_values.JxW(q_point));            // dx
              if(MANUAL_LIFTING_ON) {
                cell_rhs(i) -= eps_r * eps_0 *
                        (fe_values.shape_grad(i, q_point) *   // grad phi_i(x_q)
                          rg_gradients[q_point] *             // grad_Rg(x_q)
                          fe_values.JxW(q_point));            // dx
              }
            }
          cell->get_dof_indices(local_dof_indices);
          //for (unsigned int i = 0; i < dofs_per_cell; ++i)
          //  rhs(local_dof_indices[i]) += cell_rhs(i);
          constraints.distribute_local_to_global(cell_rhs,local_dof_indices,rhs);
        }
      }
    }

    template <int dim>
    void PrimalSolver<dim>::construct_Rg_vector() {
      AffineConstraints<double> hanging_node_constraints;
      hanging_node_constraints.clear();
      void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) =
          &DoFTools::make_hanging_node_constraints;
      // Start a side task then continue on the main thread
      Threads::Task<void> side_task =
        Threads::new_task(mhnc_p, this->dof_handler, hanging_node_constraints);

      std::map<types::global_dof_index, double> boundary_value_map;
      if (LOAD_FROM_SETUP==11) {
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 1, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 2, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 3, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 4, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
      } else {
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 1, *(this->boundary_values), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 2, *(this->boundary_values), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 9, *(this->boundary_values), boundary_value_map); //TODO: remove 9
      }
      for (const auto &boundary_value : boundary_value_map)
        this->locally_relevant_Rg_vector(boundary_value.first) = boundary_value.second;

      //cout<< std::scientific << std::setprecision(12)<<"Exact point value at EVALUATION POINT : "<<this->boundary_values->value(Point<dim>(EVALUATION_POINT_X,EVALUATION_POINT_Y))<<std::endl;

      side_task.join();
      hanging_node_constraints.close();
      hanging_node_constraints.distribute(this->locally_relevant_Rg_vector);
    }

    template <int dim>
    void PrimalSolver<dim>::interpolate_boundary_values(std::map<types::global_dof_index, double> & boundary_value_map) {
      if(MANUAL_LIFTING_ON) {
        if (LOAD_FROM_SETUP == 11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ZeroFunction<dim>(),boundary_value_map);
        } else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,Functions::ZeroFunction<dim>(),boundary_value_map);
        }
      } else {
        if (LOAD_FROM_SETUP==11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ConstantFunction<dim>(Ve),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ConstantFunction<dim>(Ve),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ConstantFunction<dim>(Vc),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ConstantFunction<dim>(Vc),boundary_value_map);
        }
        else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,*this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,*this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,*this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,*this->boundary_values,boundary_value_map);
        }
      }

    }

    // ------------------------------------------------------
    // DualSolver
    // ------------------------------------------------------

    template <int dim>
    DualSolver<dim>::DualSolver(
      parallel::distributed::Triangulation<dim> &    triangulation,
      const FiniteElement<dim> &                     fe,
      const Quadrature<dim> &                        quadrature,
      const Quadrature<dim - 1> &                    face_quadrature,
      const DualFunctional::DualFunctionalBase<dim> &dual_functional,
      const Function<dim> &                          special_rhs_function,
      const Function<dim> &                          special_boundary_values,
      const unsigned degree)
      : Base<dim>(triangulation)
      , Solver<dim>(triangulation,
                    fe,
                    quadrature,
                    face_quadrature,
                    boundary_values,
                    degree)
      , special_rhs_function(&special_rhs_function)
      , special_boundary_values(&special_boundary_values)
      , dual_functional(&dual_functional)
    {}

    template <int dim>
    void DualSolver<dim>::assemble_rhs(PETScWrappers::MPI::Vector &rhs, AffineConstraints<double> & constraints) const
    {
      dual_functional->assemble_rhs(this->dof_handler, rhs, constraints);
    }

    template <int dim>
    void DualSolver<dim>::interpolate_boundary_values(std::map<types::global_dof_index, double> & boundary_value_map) {
      if(MANUAL_LIFTING_ON) {
        if (LOAD_FROM_SETUP == 11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ZeroFunction<dim>(),boundary_value_map);
        } else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,Functions::ZeroFunction<dim>(),boundary_value_map);
        }
      } else {
        if (LOAD_FROM_SETUP==11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ZeroFunction<dim>(),boundary_value_map);
        }
        else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,this->boundary_values,boundary_value_map);
        }
      }

    }

    // Template instantiation
    template class Base<2>;
    template class Solver<2>;
    template class PrimalSolver<2>;
    template class DualSolver<2>;

  } // namespace LaplaceSolver
} // namespace IonPropulsion