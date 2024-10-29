#ifndef UTILITIES_H
#define UTILITIES_H

#include "globals.h"
#include <deal.II/base/function.h>


using namespace dealii;
constexpr double pi = 3.14159265358979323846;
std::string extract_mesh_name();

// -----------------------------------------
// RHS
// -----------------------------------------


template <int dim>
class RightHandSide0 : public Function<dim>{
public:
	virtual void value_list( 	const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
	(void)component;
	AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p)
			values[p] = 0.;
	}
};

template <int dim>
class RightHandSide1 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double k = 5.0;
		const auto x = p[0];
		const auto y = p[1];
		double piL = pi / (L*L);
		if(std::abs(x) <= k*L && y<= k*L){
			double argX = x*x * piL;
			double argY = y*y * piL;
			return eps_0 * eps_r * 2. * Ve * piL *
						 ( sin(argY) * (-2 * piL * sin(argX)*x*x + cos(argX)) +
						   sin(argX) * (-2 * piL * sin(argY)*y*y + cos(argY)) );
		}	else
			return 0.;
	}
};
template <int dim>
class RightHandSide2 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		const auto r2 = x*x + y*y;
		return 1./0.004 * eps_0 * eps_r * exp(-r2/0.004) * 4. * (1 - r2/0.004);
	}
};

template <int dim>
class RightHandSide3 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		return eps_0 * eps_r *2. / (pow(0.004,4)) * (x*x + y*y);
	}
};

template <int dim>
class RightHandSide4 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		double r2 = x*x + y*y;

    if(r2 <= Rc*Rc){
			return - eps_0 * eps_r
						 * 4. * ( AE + AD*(-2.*R*R + 4.*r2) + 3.*AC*(R*R*R*R - 4.*R*R*r2 + 3.*r2*r2));
		} else
      return 0.;
	}
};

// -----------------------------------------
// Exact solution
// -----------------------------------------


template <int dim>
class ExactSolution1 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double k = 5.0;
		const auto x = p[0];
		const auto y = p[1];
		if(std::abs(x) <= k*L && y<= k*L)
			return Ve * (1 - 
						 std::sin(std::abs(x)/L*pi) * 
						 std::sin(y/L*pi) );
						 //(1.-1./((k-1)*L)*(std::abs(x)-L)) *
						 //(1.-1./((k-1)*L)*(y-L));
		else
			return 0.;
	}
};

template <int dim>
class ExactSolution2 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		const auto r2 = x*x + y*y;
		return exp(-r2/0.004);
	}
};

template <int dim>
class ExactSolution3 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		return 1. - x*x * y*y / pow(0.004,4);
	}
};

template <int dim>
class ExactSolution4 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		
		const auto x = p[0];
		const auto y = p[1];
		double r2 = x*x + y*y;
    if(r2 <= Rc*Rc)
      return AC*(pow((r2-R*R),3)) + AD*(pow((r2-R*R),2)) + AE*(r2-R*R) + AF;
    else
      return 0.;
	}
};

#endif