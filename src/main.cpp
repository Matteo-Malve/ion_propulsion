#include <iostream>
#include "Foundamentals/Framework.h"


int main()
{
    try
    {
        // Fix dim
        const unsigned int dim = 2;

        // Define DESCRIPTOR
        Framework<dim>::ProblemDescription descriptor;
        // Fill DESCRIPTOR fields
        descriptor.primal_fe_degree = 1;
        descriptor.dual_fe_degree   = 2;
        descriptor.dual_functional = std::make_unique<EmitterFlux<dim>>();
        cout<<"\n\n\n\n\nTutto Bene\n\n\n\n\n\n";
        descriptor.max_degrees_of_freedom = 20000;
        // RUN FRAMEWORK
        Framework<dim>::run(descriptor);
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