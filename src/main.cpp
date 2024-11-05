#include "problem.h"
#include <iostream>


int main() {
	try {
		std::cout << "deal.II version: " << DEAL_II_PACKAGE_VERSION << std::endl;
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