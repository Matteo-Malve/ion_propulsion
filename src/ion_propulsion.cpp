#include "Framework.h"


int main()
{
  try
  {
    using namespace IonPropulsion;
    if (MANUAL_LIFTING_ON)
      cout<<"Lifting strategy:   Manual lifting"<<std::endl;
    else
      cout<<"Lifting strategy:   Automatic imposition of BCs"<<std::endl;

    // Describe the problem we want to solve here by passing a descriptor
    // object to the function doing the rest of the work:
    const unsigned int                 dim = 2;
    Framework<dim>::ProblemDescription descriptor;

    // First set the refinement criterion we wish to use:
    descriptor.refinement_criterion =
      Framework<dim>::ProblemDescription::dual_weighted_error_estimator;

    descriptor.primal_fe_degree = 1;
    descriptor.dual_fe_degree   = 2;

    descriptor.data = std::make_unique<Data::SetUp<Data::LogCircular_1_10<dim>, dim>>();

    const Point<dim> evaluation_point(0.0019, 0.);    // LogCircular 1:10
    //const Point<dim> evaluation_point(0.019375, 0.);    // LogCircular 1:100
    //const Point<dim> evaluation_point(0.0039, 0.0039);      // ??
    //const Point<dim> evaluation_point(0.0006, 0.0006);      // Rectangular_1_99_DeFalco PointEvaluation is sharp. Requires vertex //TODO: Extrapolation of value
    //const Point<dim> evaluation_point(0.75, 0.75);    // original-step14

    //descriptor.dual_functional = std::make_unique<DualFunctional::PointValueEvaluation<dim>>(evaluation_point);
    descriptor.dual_functional = std::make_unique<DualFunctional::StandardFluxEvaluation<dim>>(1);

    Evaluation::PointValueEvaluation<dim> postprocessor1(evaluation_point);
    Evaluation::L2_error_estimate<dim> postprocessor2(descriptor.data->get_exact_solution());
    Evaluation::H1_error_estimate<dim> postprocessor3(descriptor.data->get_exact_solution());
    Evaluation::FluxEvaluation<dim> postprocessor4;
    //Evaluation::GridOutput<dim>           postprocessor2("grid");

    descriptor.evaluator_list.push_back(&postprocessor1);
    descriptor.evaluator_list.push_back(&postprocessor2);
    descriptor.evaluator_list.push_back(&postprocessor3);
    descriptor.evaluator_list.push_back(&postprocessor4);

    // Set the maximal number of degrees of freedom after which we want the
    // program to stop refining the mesh further:
    descriptor.max_degrees_of_freedom = 2500000;

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
