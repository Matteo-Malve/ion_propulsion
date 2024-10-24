#ifndef UTILITIES_H
#define UTILITIES_H

#include "globals.h"
#include <deal.II/base/function.h>

using namespace dealii;

std::string extract_mesh_name();

template <int dim>
class RightHandSide : public Function<dim>{
public:
	virtual void value_list( 	const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
	(void)component;
	AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p)
			values[p] = 0.;
	}
};

#endif