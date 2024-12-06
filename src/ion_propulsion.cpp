#include "Framework.h"


int main()
{
  try
    {
      using namespace IonPropulsion;

      // Describe the problem we want to solve here by passing a descriptor
      // object to the function doing the rest of the work:
      const unsigned int                 dim = 2;
      Framework<dim>::ProblemDescription descriptor;

      // First set the refinement criterion we wish to use:
      descriptor.refinement_criterion =
        Framework<dim>::ProblemDescription::dual_weighted_error_estimator;
      // Here, we could as well have used <code>global_refinement</code> or
      // <code>weighted_kelly_indicator</code>. Note that the information
      // given about dual finite elements, dual functional, etc is only
      // important for the given choice of refinement criterion, and is
      // ignored otherwise.

      // Then set the polynomial degrees of primal and dual problem. We choose
      // here bi-linear and bi-quadratic ones:
      descriptor.primal_fe_degree = 1;
      descriptor.dual_fe_degree   = 2;

      // Then set the description of the test case, i.e. domain, boundary
      // values, and right hand side. These are prepackaged in classes. We
      // take here the description of <code>Exercise_2_3</code>, but you can
      // also use <code>CurvedRidges@<dim@></code>:
      descriptor.data =
        std::make_unique<Data::SetUp<Data::Exercise_2_3<dim>, dim>>();

      // Next set first a dual functional, then a list of evaluation
      // objects. We choose as default the evaluation of the value at an
      // evaluation point, represented by the classes
      // <code>PointValueEvaluation</code> in the namespaces of evaluation and
      // dual functional classes. You can also set the
      // <code>PointXDerivativeEvaluation</code> classes for the x-derivative
      // instead of the value at the evaluation point.
      //
      // Note that dual functional and evaluation objects should
      // match. However, you can give as many evaluation functionals as you
      // want, so you can have both point value and derivative evaluated after
      // each step.  One such additional evaluation is to output the grid in
      // each step.
      const Point<dim> evaluation_point(0.75, 0.75);
      descriptor.dual_functional =
        std::make_unique<DualFunctional::PointValueEvaluation<dim>>(
          evaluation_point);

      Evaluation::PointValueEvaluation<dim> postprocessor1(evaluation_point);
      Evaluation::GridOutput<dim>           postprocessor2("grid");

      descriptor.evaluator_list.push_back(&postprocessor1);
      descriptor.evaluator_list.push_back(&postprocessor2);

      // Set the maximal number of degrees of freedom after which we want the
      // program to stop refining the mesh further:
      descriptor.max_degrees_of_freedom = 20000;

      // Finally pass the descriptor object to a function that runs the entire
      // solution with it:
      Framework<dim>::run(descriptor);
    }

  // Catch exceptions to give information about things that failed:
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
