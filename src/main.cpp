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
        unsigned int grid_option = redefined_2_datafile("Load/grid_option",3);
        // Define DESCRIPTOR
        Functions::ZeroFunction <dim> zero_function;
        ProblemDescription<dim> descriptor(zero_function);
        // Fill DESCRIPTOR fields
        descriptor.primal_fe_degree = 1;
        descriptor.dual_fe_degree = 2;
        descriptor.dual_functional = std::make_unique<EmitterFlux<dim>>();
        descriptor.max_degrees_of_freedom = 20000;

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