#include <iostream>
#include "Foundamentals/PrimalSolver.h"
#include "Foundamentals/DualSolver.h"


int main()
{
    try
    {
        cout<<"--------------------------------"<<endl
            <<"Solving Primal Problem"<<endl
            <<"--------------------------------"<<endl;
        PrimalSolver<2> primal_solver;
        primal_solver.run();
        ionization_area(primal_solver.triangulation, primal_solver.dof_handler, primal_solver.solution);
        cout<<endl;

        cout<<"--------------------------------"<<endl
            <<"Solving Dual Problem"<<endl
            <<"--------------------------------"<<endl;
        DualSolver<2> dual_solver;
        dual_solver.run();
        cout<<endl;

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

