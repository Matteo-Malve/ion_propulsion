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

#ifndef LAPLACE_SOLVER_H
#define LAPLACE_SOLVER_H

#include "includes.h"
#include "Evaluation.h"
#include "DualFunctional.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{

    // ------------------------------------------------------
    // Base
    // ------------------------------------------------------
    template <int dim>
    class Base
    {
    public:
      Base(parallel::distributed::Triangulation<dim> &coarse_grid);
      virtual ~Base() = default;

      virtual void solve_problem() = 0;
      virtual void postprocess(
        const Evaluation::EvaluationBase<dim> &postprocessor) const = 0;
      virtual void         refine_grid()                            = 0;
      virtual unsigned int n_dofs() const                           = 0;

      virtual void set_refinement_cycle(const unsigned int cycle);

      virtual void output_solution()  = 0;
      virtual void update_convergence_table() = 0;
      virtual void print_convergence_table() const {};
      void print_and_reset_timer() {
        computing_timer.print_summary();
        computing_timer.reset();
      }
      TimerOutput & get_timer() {
        return computing_timer;
      }

      template <class Archive>
      void serialize(Archive &ar, const unsigned int version) {
        (void)ar;
        (void)version;
      };

      void checkpoint();
      virtual void restart() = 0;

    protected:
      const SmartPointer<parallel::distributed::Triangulation<dim>> triangulation;
      std::shared_ptr<ConvergenceTable> convergence_table;

      unsigned int refinement_cycle;

      MPI_Comm mpi_communicator;
      const unsigned int n_mpi_processes;
      const unsigned int this_mpi_process;
      ConditionalOStream pcout;

      TimerOutput        computing_timer;
    };

    // ------------------------------------------------------
    // Solver
    // ------------------------------------------------------
    template <int dim>
    class Solver : public virtual Base<dim>
    {
    public:
      Solver(parallel::distributed::Triangulation<dim> &       triangulation,
             const FiniteElement<dim> & fe,
             const Quadrature<dim> &    quadrature,
             const Quadrature<dim - 1> &face_quadrature,
             const Function<dim> &      boundary_values,
             const unsigned degree
             );
      virtual ~Solver() override;

      virtual void solve_problem() override;

      virtual void postprocess(
        const Evaluation::EvaluationBase<dim> &postprocessor) const override;
      virtual void interpolate_boundary_values(std::map<types::global_dof_index, double> &) = 0;

      virtual unsigned int n_dofs() const override;

      void update_convergence_table() override;

      void restart() override;

    protected:
      const SmartPointer<const FiniteElement<dim>>  fe;
      const SmartPointer<const Quadrature<dim>>     quadrature;
      const SmartPointer<const Quadrature<dim - 1>> face_quadrature;
      DoFHandler<dim>                               dof_handler;
      PETScWrappers::MPI::Vector                    locally_relevant_solution;
      //PETScWrappers::MPI::Vector                    homogeneous_solution;
      PETScWrappers::MPI::Vector                    locally_relevant_Rg_vector;   // TODO
      PETScWrappers::MPI::Vector                    completely_distributed_solution;
      const SmartPointer<const Function<dim>>       boundary_values;
      double                                        conservative_flux;
      MappingQ<dim>      mapping;
      IndexSet                                      locally_owned_dofs;
      IndexSet                                      locally_relevant_dofs;

      virtual void assemble_rhs(PETScWrappers::MPI::Vector &rhs, AffineConstraints<double> &) const = 0;

      virtual void construct_Rg_vector() = 0;

      virtual void retrieve_Rg() = 0;


      struct LinearSystem
      {
        LinearSystem(
          const DoFHandler<dim> &dof_handler,
          MPI_Comm & mpi_communicator,
          IndexSet & locally_owned_dofs,
          IndexSet & locally_relevant_dofs);

        void solve(
          PETScWrappers::MPI::Vector &solution,
          PETScWrappers::MPI::Vector & completely_distributed_solution) const;

        AffineConstraints<double> hanging_node_constraints;
        //SparsityPattern           sparsity_pattern;
        PETScWrappers::MPI::SparseMatrix    matrix;
        PETScWrappers::MPI::Vector          rhs;
        PETScWrappers::MPI::SparseMatrix    Umatrix;
      };

      std::unique_ptr<LinearSystem>                 linear_system_ptr;

    private:

      struct AssemblyScratchData
      {
        AssemblyScratchData(const FiniteElement<dim> &fe,
                            const Quadrature<dim> &   quadrature,
                            MappingQ<dim> & mapping);
        AssemblyScratchData(const AssemblyScratchData &scratch_data);

        MappingQ<dim>  mapping;
        FEValues<dim> fe_values;
      };

      struct AssemblyCopyData
      {
        FullMatrix<double>                   cell_matrix;
        std::vector<types::global_dof_index> local_dof_indices;
      };


      void assemble_linear_system(LinearSystem &linear_system);

      void local_assemble_matrix(
        const typename DoFHandler<dim>::active_cell_iterator &cell,
        AssemblyScratchData &                                 scratch_data,
        AssemblyCopyData &                                    copy_data) const;


      void copy_local_to_global(const AssemblyCopyData &copy_data,
                                LinearSystem &          linear_system) const;


    };


    // ------------------------------------------------------
    // PrimalSolver
    // ------------------------------------------------------
    template <int dim>
    class PrimalSolver : public Solver<dim>
    {
    public:
      PrimalSolver(parallel::distributed::Triangulation<dim> &       triangulation,
                   const FiniteElement<dim> & fe,
                   const Quadrature<dim> &    quadrature,
                   const Quadrature<dim - 1> &face_quadrature,
                   const Function<dim> &      rhs_function,
                   const Function<dim> &      boundary_values,
                   const unsigned degree);

      virtual void output_solution()  override;

    protected:
      const SmartPointer<const Function<dim>>       rhs_function;
      double                                        conservative_flux;
      Vector<double>                                Au;

      virtual void assemble_rhs(PETScWrappers::MPI::Vector &rhs
                                      , AffineConstraints<double> &) const override;
      virtual void construct_Rg_vector() override;
      void interpolate_boundary_values(std::map<types::global_dof_index, double> &) override;

    private:
      void retrieve_Rg() override {
        this->locally_relevant_solution += this->locally_relevant_Rg_vector;
      }
    };

    // ------------------------------------------------------
    // DualSolver
    // ------------------------------------------------------
    template <int dim>
    class DualSolver : public Solver<dim>
    {
    public:
      DualSolver(
        parallel::distributed::Triangulation<dim> &    triangulation,
        const FiniteElement<dim> &                     fe,
        const Quadrature<dim> &                        quadrature,
        const Quadrature<dim - 1> &                    face_quadrature,
        const DualFunctional::DualFunctionalBase<dim> &dual_functional,
        const Function<dim> &                          special_rhs_function,
        const Function<dim> &                          special_boundary_values,
        const unsigned degree);

    protected:
      const SmartPointer<const Function<dim>>       special_rhs_function;
      const SmartPointer<const Function<dim>>       special_boundary_values;

      const SmartPointer<const DualFunctional::DualFunctionalBase<dim>>
                   dual_functional;
      virtual void assemble_rhs(PETScWrappers::MPI::Vector &rhs
                                      , AffineConstraints<double> &) const override;
      static const Functions::ZeroFunction<dim> boundary_values;

      /*void assemble_conservative_flux_rhs(
        Vector<double> &rhs
        ,SparseMatrix<double> & Umatrix
        );*/

    private:
      virtual void construct_Rg_vector() override {};
      void retrieve_Rg() override {};
      void interpolate_boundary_values(std::map<types::global_dof_index, double> &) override;

    };

    template <int dim>
    const Functions::ZeroFunction<dim> DualSolver<dim>::boundary_values;



  } // namespace LaplaceSolver
}

#endif