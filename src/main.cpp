#include "problem.h"
#include <iostream>


int main() {
	try {
		std::cout << std::endl;
		std::cout << "deal.II version: " << DEAL_II_PACKAGE_VERSION << std::endl;

		std::cout << "Loaded Mesh name: " << extract_mesh_name() << std::endl;

		std::cout << "Chosen goal functional: " << GOAL_FUNCTIONAL << std::endl;

		std::string onoff = ENABLE_CONVERGENCE_ANALYSIS ? "ON" : "OFF";
		std::cout << "Convergence analysis: " << onoff << std::endl;

		if(GOAL_FUNCTIONAL=="PointValue" || GOAL_FUNCTIONAL=="PointYDerivative")
			cout<<"Evalaution point = ("<<EVALUATION_POINT(0)<<" , "<<EVALUATION_POINT(1)<<")"<<endl;

		Problem<2> iprop_problem;
		iprop_problem.run();
	} catch (std::exception &exc) {
		std::cerr << "\nException on processing: " << exc.what() << "\nAborting!\n";
		return 1;
	} catch (...) {
		std::cerr << "\nUnknown exception!\nAborting!\n";
		return 1;
	}
	return 0;
}