#include <iostream>
#include "../include/Problem.h"
//#include "Problem.h"
int main()
{
    try
    {
        // Il numero in <> indica la dimensione del problema

        Problem<2> wire_poisson_2d; // Definisco il tipo di problema...
        wire_poisson_2d.run(); // ... e lo simulo

        // Non preoccupatevi del caso 3D: non funziona...
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

        return 1; // Report an error
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

