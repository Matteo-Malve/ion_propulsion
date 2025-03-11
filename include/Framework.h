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

#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#include "Evaluation.h"
#include "Data.h"
#include "DualFunctional.h"
#include "Refinement.h"

namespace IonPropulsion{
	using namespace dealii;

	template <int dim>
  struct Framework
  {
  public:
    using Evaluator     = Evaluation::EvaluationBase<dim>;
    using EvaluatorList = std::list<Evaluator *>;

    struct ProblemDescription
    {
      unsigned int primal_fe_degree;
      unsigned int dual_fe_degree;
      std::unique_ptr<const Data::SetUpBase<dim>> data;

      enum RefinementCriterion
      {
        dual_weighted_error_estimator,
        global_refinement,
        kelly_indicator,
        weighted_kelly_indicator
      };

      RefinementCriterion refinement_criterion;

      std::unique_ptr<const DualFunctional::DualFunctionalBase<dim>>
        dual_functional;

      EvaluatorList evaluator_list;

      std::unique_ptr<const Function<dim>> kelly_weight;

      unsigned int max_degrees_of_freedom;

      unsigned int mapping_degree;

      ProblemDescription();
    };

    static void run(const ProblemDescription &descriptor, const bool do_restart);
  };

}


#endif