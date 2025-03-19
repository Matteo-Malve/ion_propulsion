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

#include "Framework.h"

namespace IonPropulsion{
  using namespace dealii;

  template <int dim>
  Framework<dim>::ProblemDescription::ProblemDescription()
    : primal_fe_degree(1)
    , dual_fe_degree(2)
    , refinement_criterion(dual_weighted_error_estimator)
    , max_degrees_of_freedom(20000)
    , mapping_degree(2)
  {}

  template <int dim>
  void Framework<dim>::run(const ProblemDescription &descriptor)
  {
    MPI_Comm mpi_communicator(MPI_COMM_WORLD);
    //const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator));
    const unsigned int this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator));
    ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

    parallel::distributed::Triangulation<dim> triangulation(
      mpi_communicator,
      typename Triangulation<dim>::MeshSmoothing(
        Triangulation<dim>::smoothing_on_refinement |
        Triangulation<dim>::smoothing_on_coarsening)
        );
    descriptor.data->create_coarse_grid(triangulation);

    // then a set of finite elements and appropriate quadrature formula:
    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);

    // Next, select one of the classes implementing different refinement
    // criteria.
    std::unique_ptr<LaplaceSolver::Base<dim>> solver;
    switch (descriptor.refinement_criterion)
      {
        case ProblemDescription::dual_weighted_error_estimator:
        {
          solver = std::make_unique<LaplaceSolver::WeightedResidual<dim>>(
            triangulation,
            primal_fe,
            dual_fe,
            quadrature,
            face_quadrature,
            descriptor.data->get_right_hand_side(),
            descriptor.data->get_boundary_values(),
            *descriptor.dual_functional,
            descriptor.mapping_degree);
          break;
          }
        case ProblemDescription::global_refinement:
          {
            solver = std::make_unique<LaplaceSolver::RefinementGlobal<dim>>(
              triangulation,
              primal_fe,
              quadrature,
              face_quadrature,
              descriptor.data->get_right_hand_side(),
              descriptor.data->get_boundary_values(),
              descriptor.mapping_degree);
            break;
          }
        default:
          AssertThrow(false, ExcInternalError());
      }

    for (unsigned int step = 0; true; ++step)
      {

      // ------------------------------------------------------
      // Solution
      // ------------------------------------------------------

      pcout << "Refinement cycle: " << step << std::endl;
      solver->set_refinement_cycle(step);

      solver->solve_problem();

      MPI_Barrier(mpi_communicator);
      solver->update_convergence_table();

      MPI_Barrier(mpi_communicator);
      solver->output_solution();

      MPI_Barrier(mpi_communicator);
      pcout << "   Number of degrees of freedom=" << solver->n_dofs() << std::endl;


      // ------------------------------------------------------
      // Evaluation
      // ------------------------------------------------------
      {
        TimerOutput::Scope t(solver->get_timer(), "evaluation");
        for (const auto &evaluator : descriptor.evaluator_list)
        {
          evaluator->set_refinement_cycle(step);
          solver->postprocess(*evaluator);
        }
      }

      // ------------------------------------------------------
      // Refinement
      // ------------------------------------------------------
      MPI_Barrier(mpi_communicator);
      unsigned int DoFs_before_refinement = solver->n_dofs();
      solver->refine_grid();

      // ------------------------------------------------------
      // Output
      // ------------------------------------------------------
      MPI_Barrier(mpi_communicator);
      if (this_mpi_process == 0) {
        solver->print_convergence_table();
        CSVLogger::getInstance().flushRow();
      }
      MPI_Barrier(mpi_communicator);
      //if (REFINEMENT_CRITERION==2)
        //solver->checkpoint();

      // ------------------------------------------------------
      // Closure
      // ------------------------------------------------------
      MPI_Barrier(mpi_communicator);
      cout<<std::defaultfloat;
      MPI_Barrier(mpi_communicator);
      solver->print_and_reset_timer();
      MPI_Barrier(mpi_communicator);

      if (DoFs_before_refinement > descriptor.max_degrees_of_freedom)
        break;

    }

  // Clean up the screen after the loop has run:
  std::cout << std::endl;
}

  // Template instantiation
  template struct Framework<2>;
} // namespace IonPropulsion