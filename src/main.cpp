#include <iostream>
#include "Foundamentals/Framework.h"

static GetPot redefined_2_datafile("../data_setup");


int main()
{
    try
    {
        // Fix dim
        const unsigned int dim = 2;
        // Fix grid_option
        unsigned int grid_option = redefined_2_datafile("grid_option",3 );
        // Define DESCRIPTOR
        Functions::ZeroFunction<dim> zero_function;
        ProblemDescription<dim> descriptor(zero_function);
        // Fill DESCRIPTOR fields
        descriptor.primal_fe_degree = redefined_2_datafile("primal_fe_degree",5);
        descriptor.dual_fe_degree = redefined_2_datafile("dual_fe_degree",5);
        descriptor.dual_functional = std::make_unique<EmitterFlux<dim>>();
        descriptor.max_degrees_of_freedom = redefined_2_datafile("max_degrees_of_freedom",3000);
        descriptor.max_number_of_refinements = redefined_2_datafile("max_number_of_refinements",5);
        cout<< "[Main] Data setup:" <<endl
            << "        - primal_fe_degree                   "<<descriptor.primal_fe_degree<<endl
            << "        - dual_fe_degree                     "<<descriptor.dual_fe_degree<<endl
            << "        - max_degrees_of_freedom             "<<descriptor.max_degrees_of_freedom<<endl
            << "        - max_number_of_refinements          "<<descriptor.max_number_of_refinements<<endl;

        // RUN FRAMEWORK
        framework_run(descriptor,grid_option);

    }

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