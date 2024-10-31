#ifndef UTILITIES_H
#define UTILITIES_H

#include "globals.h"
#include <deal.II/base/function.h>


using namespace dealii;
constexpr double pi = 3.14159265358979323846;
std::string extract_mesh_name();

// -----------------------------------------
// 0
// -----------------------------------------


template <int dim>
class RightHandSide0 : public Function<dim>{
public:
	virtual void value_list( 	const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
	(void)component;
	AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			const auto x = point_list[p][0];
			const auto y = point_list[p][1];
			double r = sqrt(x*x + y*y);
			values[p] = -Ve;
		}
	}
};

// -----------------------------------------
// 1 - broken
// -----------------------------------------

template <int dim>
class RightHandSide1 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double k = 5.0;
		const auto x = p[0];
		const auto y = p[1];
		double piL = pi / (l*l);
		if(std::abs(x) <= k*l && y<= k*l){
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
class ExactSolution1 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double k = 5.0;
		const auto x = p[0];
		const auto y = p[1];
		if(std::abs(x) <= k*l && y<= k*l)
			return Ve * (1 - 
						 std::sin(std::abs(x)/l*pi) * 
						 std::sin(y/l*pi) );
						 //(1.-1./((k-1)*l)*(std::abs(x)-l)) *
						 //(1.-1./((k-1)*l)*(y-l));
		else
			return 0.;
	}
};

// -----------------------------------------
// 2 - broken
// -----------------------------------------

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

// -----------------------------------------
// 3 - broken
// -----------------------------------------

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
class ExactSolution3 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		return 1. - x*x * y*y / pow(0.004,4);
	}
};

// -----------------------------------------
// 4 - WORKS
// -----------------------------------------

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
// 5 - to test - for rectangular
// -----------------------------------------


template <int dim>
class ExactSolution5 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		
		double sigma2 = 0.0000005;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
    
		return exp(-(r-R)*(r-R)/(2.*sigma2))
					 * Ve
					 * (1. - sin(x/l*pi) * sin(y/l*pi));
      
	}
};

template <int dim>
class RightHandSide5 : public Function<dim>{
public:
	virtual void value_list( 	const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		cout<<"In value_list"<<endl;
		AssertDimension (point_list.size(), values.size()); // Size check

		double sigma2 = 0.0000005;
		for (unsigned int p=0; p<point_list.size(); ++p){
			const auto x = point_list[p][0];
			const auto y = point_list[p][1];
			double r = sqrt(x*x + y*y);
			double argX = pi * x / l;
			double argY = pi * y / l;

			// Break down the formula for readability
			double expTerm = std::exp(-std::pow(R - r, 2) / (2 * sigma2));
		
			double sinArgX = std::sin(argX);
			double sinArgY = std::sin(argY);
			double cosArgX = std::cos(argX);
			double cosArgY = std::cos(argY);
			
			double term1 = l * l * R * R * r;
			double term2 = -l * l * R * R * r * sinArgX * sinArgY;
			double term3 = -l * l * R * sigma2 * sinArgX * sinArgY;
			double term4 = l * l * R * sigma2;
			double term5 = 2 * l * l * R * x * x * sinArgX * sinArgY;
			double term6 = -2 * l * l * R * x * x;
			double term7 = 2 * l * l * R * y * y * sinArgX * sinArgY;
			double term8 = -2 * l * l * R * y * y;
			double term9 = -2 * l * l * sigma2 * r;
			double term10 = 2 * l * l * sigma2 * r * sinArgX * sinArgY;
			double term11 = l * l * x * x * r;
			double term12 = l * l * y * y * r;
			double term13 = -l * l * x * x * r * sinArgX * sinArgY;
			double term14 = -l * l * y * y * r * sinArgX * sinArgY;
			double term15 = 2 * pi * l * sigma2 * y * sinArgX * cosArgY * (r - R);
			double term16 = 2 * pi * l * sigma2 * x * cosArgX * sinArgY * (r - R);
			double term17 = 2 * pi * pi * sigma2 * sigma2 * r * sinArgX * sinArgY;
			
			double numerator = expTerm * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 +
																		term9 + term10 + term11 + term12 + term13 + term14 + term15 +
																		term16 + term17);
			double denominator = l * l * sigma2 * sigma2 * r;
			double laplacian = numerator / denominator;
			values[p] = - eps_0 * eps_r * laplacian;
		}
			
		
	}
};

// -----------------------------------------
// 6 - to test - paraola per cerchi concentrici
// -----------------------------------------


template <int dim>
class RightHandSide6 : public Function<dim>{
public:
	virtual void value_list( 	const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
	(void)component;
	AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			const auto x = point_list[p][0];
			const auto y = point_list[p][1];
			double r = sqrt(x*x + y*y);
			values[p] = - eps_0 * eps_r *
									Ve * ( l + L - 4.*r ) / r;
		}
	}
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);

		return - eps_0 * eps_r *
					 Ve * ( l + L - 4.*r ) / r;
	}
};

template <int dim>
class ExactSolution6 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);

		return - Ve * (r-l) * (r-L);
	}
};

#endif