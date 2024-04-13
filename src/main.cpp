#include <iostream>
#include "Problem.h"

int main()
{
    try
    {
        Constants constants;
        constants.eps_0 = 8.854e-12; // [F/m]
        constants.eps_r = 1.0006;
        constants.Vmax = 2.0e+4; // [V]
        constants.E_ON = 3.31e+6; // [V/m]
        constants.R = 2.5e-4; // [m]
        constants.nn = 2;
        constants.L = constants.nn * constants.R; // [m]
        constants.X = 0.0; // [m]
        constants.g = 0.2; // [m]
        constants.mesh_height = 0.1; // [m]

	    Problem<2> problem(constants);
	    problem.run();
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
