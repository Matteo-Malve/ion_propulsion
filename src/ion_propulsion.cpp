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
      descriptor.refinement_criterion = Framework<dim>::ProblemDescription::dual_weighted_error_estimator;
      //descriptor.refinement_criterion = Framework<dim>::ProblemDescription::global_refinement;

      descriptor.primal_fe_degree = 1;
      descriptor.dual_fe_degree   = 2;

      descriptor.data = std::make_unique<Data::SetUp<Data::Exercise_2_3<dim>, dim>>();
      //descriptor.data = std::make_unique<Data::SetUp<Data::Rectangle_1_99<dim>, dim>>();

      const Point<dim> evaluation_point(0.75, 0.75);
      //const Point<dim> evaluation_point(0.0039, 0.0039);
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
