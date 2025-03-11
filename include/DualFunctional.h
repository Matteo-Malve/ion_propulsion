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

#ifndef DUALFUNCTIONAL_H
#define DUALFUNCTIONAL_H

#include "includes.h"


namespace IonPropulsion{

	using namespace dealii;

  namespace DualFunctional
  {
    // ------------------------------------------------------
    // DualFunctionalBase
    // ------------------------------------------------------
    template <int dim>
    class DualFunctionalBase : public Subscriptor
    {
    protected:
      MappingQ<dim>      mapping;
    public:
      DualFunctionalBase(const unsigned mapping_degree);
      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const = 0;
      /*virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs,
                                std::unique_ptr<LaplaceSolver::Base<dim>> &) const {};*/
    };
    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------
    template <int dim>
    class PointValueEvaluation : public DualFunctionalBase<dim>
    {
    public:
      PointValueEvaluation(
        const unsigned mapping_degree,
        const Point<dim> &evaluation_point);

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    protected:
      const Point<dim> evaluation_point;
    };

    // ------------------------------------------------------
    // PointXDerivativeEvaluation
    // ------------------------------------------------------
    template <int dim>
    class PointXDerivativeEvaluation : public DualFunctionalBase<dim>
    {
    public:
      PointXDerivativeEvaluation(
        const unsigned mapping_degree,
        const Point<dim> &evaluation_point);

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    protected:
      const Point<dim> evaluation_point;
    };

    // ------------------------------------------------------
    // StandardFluxEvaluation
    // ------------------------------------------------------
    template <int dim>
    class StandardFluxEvaluation : public DualFunctionalBase<dim>
    {
    public:
      StandardFluxEvaluation(
        const unsigned mapping_degree,
        const std::set<unsigned int> &boundary_ids);

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const override;

    protected:
      const std::set<unsigned int> boundary_ids;
    };



  } // namespace DualFunctional

}
#endif