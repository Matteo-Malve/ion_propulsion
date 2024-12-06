#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#include "Evaluation.h"
#include "Data.h"
#include "DualFunctional.h"
#include "LaplaceSolver.h"

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

      ProblemDescription();
    };

    static void run(const ProblemDescription &descriptor);
  };

}


#endif