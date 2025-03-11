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

#ifndef EVALUATION_H
#define EVALUATION_H

#include "includes.h"
#include "CSVLogger.h"


namespace IonPropulsion{
  using namespace dealii;
  namespace Evaluation{

    // ------------------------------------------------------
    // EvaluationBase
    // ------------------------------------------------------

    template <int dim>
    class EvaluationBase
    {
    public:
      EvaluationBase(const unsigned degree);
      virtual ~EvaluationBase() = default;

      void set_refinement_cycle(const unsigned int refinement_cycle);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const PETScWrappers::MPI::Vector & solution,
                                                        const parallel::distributed::Triangulation<dim> &       triangulation) const = 0;

    protected:
      unsigned int refinement_cycle;
      MappingQ<dim>      mapping;
      MPI_Comm mpi_communicator;
      const unsigned int n_mpi_processes;
      const unsigned int this_mpi_process;
      ConditionalOStream pcout;

    };

    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------

    template <int dim>
    class PointValueEvaluation : public EvaluationBase<dim>
    {
    public:
      PointValueEvaluation(const unsigned degree, const Point<dim> &evaluation_point);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const PETScWrappers::MPI::Vector & solution,
                                                        const parallel::distributed::Triangulation<dim> &       triangulation) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    private:
      const Point<dim> evaluation_point;
    };

    // ------------------------------------------------------
    // FluxEvaluation
    // ------------------------------------------------------

    template <int dim>
    class FluxEvaluation : public EvaluationBase<dim>
    {
    public:
      FluxEvaluation(const unsigned degree, const std::set<unsigned int> &boundary_ids);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const PETScWrappers::MPI::Vector & solution,
                                                        const parallel::distributed::Triangulation<dim> &       triangulation) const override;

    private:
      const std::set<unsigned int> boundary_ids;

    };


    // ------------------------------------------------------
    // GridOutput
    // ------------------------------------------------------

    template <int dim>
    class GridOutput : public EvaluationBase<dim>
    {
    public:
      GridOutput(const unsigned degree, const std::string &output_name_base);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const PETScWrappers::MPI::Vector & solution,
                                                        const parallel::distributed::Triangulation<dim> &       triangulation) const override;

    private:
      const std::string output_name_base;

    };

    // ------------------------------------------------------
    // L2_error_estimate
    // ------------------------------------------------------

    template <int dim>
    class L2_error_estimate : public EvaluationBase<dim>
    {
    public:
      L2_error_estimate(const unsigned degree, const Function<dim> & analytical_solution);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const PETScWrappers::MPI::Vector & solution,
                                                        const parallel::distributed::Triangulation<dim> &       triangulation) const override;

    private:
      const std::string output_name_base;
      const SmartPointer<const Function<dim>> analytical_solution;
    };

    // ------------------------------------------------------
    // H1_error_estimate
    // ------------------------------------------------------

    template <int dim>
    class H1_error_estimate : public EvaluationBase<dim>
    {
    public:
      H1_error_estimate(const unsigned degree, const Function<dim> & analytical_solution);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const PETScWrappers::MPI::Vector & solution,
                                                        const parallel::distributed::Triangulation<dim> &       triangulation) const override;

    private:
      const std::string output_name_base;
      const SmartPointer<const Function<dim>> analytical_solution;
    };

  } // namespace Evaluation
}

#endif