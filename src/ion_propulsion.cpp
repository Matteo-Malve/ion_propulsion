#include "Framework.h"


int main(int argc, char **argv)
{
  try
  {
    using namespace IonPropulsion;

    const std::string configFile = (argc > 1) ? argv[1] : "../config.yaml";
    GlobalConstants::initialize(configFile);
    std::cout << "Output path: " << OUTPUT_PATH << std::endl;

    useGlobalConstants();
    printParsedConstants();

    // Describe the problem we want to solve here by passing a descriptor
    // object to the function doing the rest of the work:
    const unsigned int                 dim = 2;
    Framework<dim>::ProblemDescription descriptor;

    // First set the refinement criterion we wish to use:
    if (REFINEMENT_CRITERION==1)
      descriptor.refinement_criterion = Framework<dim>::ProblemDescription::global_refinement;
    else if (REFINEMENT_CRITERION==2 || REFINEMENT_CRITERION==3)
       descriptor.refinement_criterion = Framework<dim>::ProblemDescription::dual_weighted_error_estimator;
    else
      DEAL_II_NOT_IMPLEMENTED();

    descriptor.primal_fe_degree = 1;
    descriptor.dual_fe_degree   = 2;

    descriptor.mapping_degree = MAPPING_DEGREE;

    if(LOAD_FROM_SETUP == 0)
      descriptor.data = std::make_unique<Data::SetUp<Data::SetupNone<dim>, dim>>();
    else if(LOAD_FROM_SETUP == 1)
      descriptor.data = std::make_unique<Data::SetUp<Data::CurvedRidges<dim>, dim>>();
    else if(LOAD_FROM_SETUP == 2)
      descriptor.data = std::make_unique<Data::SetUp<Data::Exercise_2_3<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 3)
      descriptor.data = std::make_unique<Data::SetUp<Data::Rectangle_1_99<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 4)
      descriptor.data = std::make_unique<Data::SetUp<Data::LogCircular_1_10<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 5)
      descriptor.data = std::make_unique<Data::SetUp<Data::LogCircular_1_100<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 6)
      descriptor.data = std::make_unique<Data::SetUp<Data::Rectangle_1_99_manifold<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 7)
      descriptor.data = std::make_unique<Data::SetUp<Data::angle_step14_forced<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 8)
      descriptor.data = std::make_unique<Data::SetUp<Data::angle_Rectangle_1_100_forced<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 9)
      descriptor.data = std::make_unique<Data::SetUp<Data::CircularStep14<dim>, dim>>();
    else if (LOAD_FROM_SETUP == 10)
      descriptor.data = std::make_unique<Data::SetUp<Data::LogCircular_1_2<dim>, dim>>();
    else
      DEAL_II_NOT_IMPLEMENTED();

    const Point<dim> evaluation_point(EVALUATION_POINT_X,EVALUATION_POINT_Y);

    /*
    //const Point<dim> evaluation_point(0.0019, 0.);    // LogCircular 1:10
    const Point<dim> evaluation_point(0.019375, 0.);    // LogCircular 1:100
    //const Point<dim> evaluation_point(0.0039, 0.0039);      // ??
    //const Point<dim> evaluation_point(0.0006, 0.0006);      // Rectangular_1_99_DeFalco PointEvaluation is sharp. Requires vertex //TODO: Extrapolation of value
    //const Point<dim> evaluation_point(0.75, 0.75);    // original-step14
    */

    if (REFINEMENT_CRITERION>1) {
      if (DUAL_FUNCTIONAL==1)
        descriptor.dual_functional = std::make_unique<DualFunctional::PointValueEvaluation<dim>>(MAPPING_DEGREE,evaluation_point);
      else if (DUAL_FUNCTIONAL==2 || DUAL_FUNCTIONAL==3)
        descriptor.dual_functional = std::make_unique<DualFunctional::StandardFluxEvaluation<dim>>(MAPPING_DEGREE,1);
      else
        DEAL_II_NOT_IMPLEMENTED();
    }


    Evaluation::PointValueEvaluation<dim> postprocessor1(descriptor.mapping_degree,evaluation_point);
    Evaluation::L2_error_estimate<dim> postprocessor2(descriptor.mapping_degree,descriptor.data->get_exact_solution());
    Evaluation::H1_error_estimate<dim> postprocessor3(descriptor.mapping_degree,descriptor.data->get_exact_solution());
    Evaluation::FluxEvaluation<dim> postprocessor4(descriptor.mapping_degree);
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
