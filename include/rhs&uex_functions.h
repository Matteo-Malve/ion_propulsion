#ifndef RHS_UEX_FUNCTIONS_H
#define RHS_UEX_FUNCTIONS_H

#include "globals.h"
#include <deal.II/base/function.h>


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
			values[p] = -Ve;
		}
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
      return AC*(pow((r2-Ri*Ri),3)) + AD*(pow((r2-Ri*Ri),2)) + AE*(r2-Ri*Ri) + AF;
    else
      return 0.;
	}

	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	}

	virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
		(void)component;
		
		const auto x = p[0];
		const auto y = p[1];
		double r2 = x*x + y*y;

		Tensor<1, dim> grad;

		grad[0] = (r2 <= Rc*Rc) ?   2*x * ( AE + 2*AD * (r2 - Ri*Ri) + 3*AC * (r2 - Ri*Ri)*(r2 - Ri*Ri))    : 0.;
		grad[1] = (r2 <= Rc*Rc) ?   2*y * ( AE + 2*AD * (r2 - Ri*Ri) + 3*AC * (r2 - Ri*Ri)*(r2 - Ri*Ri))    : 0.;
		return grad;
	};

	double emitter_flux() const{
		if(PATH_TO_MESH == "../mesh/cerchi_concentrici.msh")
			return - 4 * pi * Ri*Ri * AE;
		else{
			cout<<"!! No manual exact computation of the emitter flux has been provided for this mesh."<<endl;
			abort();
		}
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
						 * 4. * ( AE + AD*(-2.*Ri*Ri + 4.*r2) + 3.*AC*(Ri*Ri*Ri*Ri - 4.*Ri*Ri*r2 + 3.*r2*r2));
		} else
      return 0.;
	}

	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
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

	virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
		(void)component;
		double sigma2 = 0.0000005;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double argX = pi * x / l;
		double argY = pi * y / l;

		double expTerm = std::exp(- std::pow(R - r, 2) / (2 * sigma2));
	
		double sinArgX = std::sin(argX);
		double sinArgY = std::sin(argY);
		double cosArgX = std::cos(argX);
		double cosArgY = std::cos(argY);


		Tensor<1, dim> grad;

		grad[0] = - expTerm * Ve * ( -pi * sigma2 * r * cosArgX * sinArgY   +    L * x * (r-R) * (sinArgX*sinArgY - 1)) 
							/ (L * sigma2 * r);
		grad[1] = expTerm * Ve * ( -pi * sigma2 * r * cosArgY * sinArgX   +    L * y * (r-R) * (sinArgX*sinArgY - 1)) 
							/ (L * sigma2 * r);
		return grad;
	}	

};

template <int dim>
class GradX : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double sigma2 = 0.0000005;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double argX = pi * x / l;
		double argY = pi * y / l;

		double expTerm = std::exp(- std::pow(R - r, 2) / (2 * sigma2));
	
		double sinArgX = std::sin(argX);
		double sinArgY = std::sin(argY);
		double cosArgX = std::cos(argX);

		return -expTerm * Ve * ( -pi * sigma2 * r * cosArgX * sinArgY   +    L * x * (r-R) * (sinArgX*sinArgY - 1)) 
							/ (L * sigma2 * r);
	}
};
template <int dim>
class GradY : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double sigma2 = 0.0000005;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double argX = pi * x / l;
		double argY = pi * y / l;

		double expTerm = std::exp(- std::pow(R - r, 2) / (2 * sigma2));
	
		double sinArgX = std::sin(argX);
		double sinArgY = std::sin(argY);
		double cosArgY = std::cos(argY);
		
		return expTerm * Ve * ( -pi * sigma2 * r * cosArgY * sinArgX   +    L * y * (r-R) * (sinArgX*sinArgY - 1)) 
							/ (L * sigma2 * r);
	}
};

template <int dim>
class RightHandSide5 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;

		double sigma2 = 0.0000005;
		
		const auto x = p[0];
		const auto y = p[1];
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
		double laplacian = Ve * numerator / denominator;
		return - eps_0 * eps_r * laplacian;
	}
	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	}
};

// -----------------------------------------
// 5b 
// -----------------------------------------


template <int dim>
class ExactSolution5b : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		
		double sigma2 = 0.0000005;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
		double freq = 0.5;
		double argX = (std::abs(x)-l) * pi / l * freq;
		double argY = (std::abs(y)-l) * pi / l * freq;
    
		return expTerm * Ve * (1 - sin(argX) * sin(argY));
      
	}

	virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
		(void)component;
		double sigma2 = 0.0000005;
		double freq = 0.5;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
		double argX = (std::abs(x)-l) * pi / l * freq;
		double argY = (std::abs(y)-l) * pi / l * freq;
		double signX = x>=0 ? +1. : -1.;
    double signY = y>=0 ? +1. : -1.;

		Tensor<1, dim> grad;

		grad[0] = Ve * expTerm * ( (- 1/sigma2 ) * (r-l) * (x/r) * (1 - sin(argX) * sin(argY))
				 												-  cos(argX) * sin(argY) * (pi/l*freq * signX));
		grad[1] = Ve * expTerm * ( (- 1/sigma2 ) * (r-l) * (y/r) * (1 - sin(argX) * sin(argY))
				 												- sin(argX) * cos(argY) * (pi/l*freq * signY));
		return grad;
};
};

template <int dim>
class GradXb : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double sigma2 = 0.0000005;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
		double freq = 0.5;
		double argX = (std::abs(x)-l) * pi / l * freq;
		double argY = (std::abs(y)-l) * pi / l * freq;
		double signX = x>=0 ? +1. : -1.;
    
		return Ve * expTerm * (- 1/sigma2 ) * (r-l) * (x/r) * (1 - sin(argX) * sin(argY))
				 - Ve * expTerm * cos(argX) * sin(argY) * (pi/l*freq * signX);
	}
};
template <int dim>
class GradYb : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double sigma2 = 0.0000005;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
		double freq = 0.5;
		double argX = (std::abs(x)-l) * pi / l * freq;
		double argY = (std::abs(y)-l) * pi / l * freq;
    double signY = y>=0 ? +1. : -1.;
    
		return Ve * expTerm * (- 1/sigma2 ) * (r-l) * (y/r) * (1 - sin(argX) * sin(argY))
				 - Ve * expTerm * sin(argX) * cos(argY) * (pi/l*freq * signY);
	}
};


template <int dim>
class RightHandSide5b : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		double sigma2 = 0.0000005;
		double freq = 0.5;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);
		double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
		double argX = (std::abs(x)-l) * pi / l * freq;
		double argY = (std::abs(y)-l) * pi / l * freq;
		double signX = x>=0 ? +1. : -1.;
    double signY = y>=0 ? +1. : -1.;

		return - eps_0 * eps_r * 
						(1 / ( l*l * sigma2*sigma2 * r)) * expTerm * Ve *
						(   2 * freq * l * pi * sigma2 * x * (r-l) * cos(argX) * sin(argY) * signX +
								freq*freq * pi*pi * sigma2*sigma2 * r * sin(argX) * sin(argY) +
								2 * freq * l * pi * sigma2 * y * (r-l) * cos(argY) * sin(argX) * signY +
								freq*freq * pi*pi * sigma2*sigma2 * r * sin(argX) * sin(argY) -
								l * ( l * ( l*l * r + r * (r*r - 2*sigma2) +l * (sigma2 - 2*r*r) ) *  (- 1 + sin(argX) * sin(argY) ))
						);
	};
	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	};
};

// -----------------------------------------
// 6 - to test - parabola per cerchi concentrici
// -----------------------------------------


template <int dim>
class RightHandSide6 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);

		return - eps_0 * eps_r *
					 Ve * ( l + L - 4.*r ) / r;
	}

	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	};
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
	};
	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			const auto x = point_list[p][0];
			const auto y = point_list[p][1];
			double r = sqrt(x*x + y*y);

			values[p] = - Ve * (r-l) * (r-L);
		};
	};

	virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);

		Tensor<1, dim> grad;

		grad[0] = 2 * Ve * x * ( l/r - 1);
		grad[1] = 2 * Ve * y * ( l/r - 1);
		return grad;
};
};

// -----------------------------------------
// 7 - Gauss * Gauss
// -----------------------------------------


template <int dim>
class RightHandSide7 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;

		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);

		return - eps_0 * eps_r *
					 Ve * ( l + L - 4.*r ) / r;
	}

	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	};
};

template <int dim>
class ExactSolution7 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);

		return - Ve * (r-l) * (r-L);
	};
	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			const auto x = point_list[p][0];
			const auto y = point_list[p][1];
			double r = sqrt(x*x + y*y);

			values[p] = - Ve * (r-l) * (r-L);
		};
	};

	virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		double r = sqrt(x*x + y*y);

		Tensor<1, dim> grad;

		grad[0] = 2 * Ve * x * ( l/r - 1);
		grad[1] = 2 * Ve * y * ( l/r - 1);
		return grad;
};
};

// -----------------------------------------
// 8 - sinsin
// -----------------------------------------

template <int dim>
class ExactSolution8 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];

		return sin(pi/L*x)*sin(pi/L*y);
	};
	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	};

	virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		

		Tensor<1, dim> grad;

		grad[0] = pi * cos(pi/L*x) * sin(pi/L*y) / L;
		grad[1] = pi * sin(pi/L*x) * cos(pi/L*y) / L;

		return grad;
};
};


template <int dim>
class RightHandSide8 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];

		double laplacian = - pi*pi / (L*L) * sin(pi/L*x) * sin(pi/L*y);
		return - eps_0 * eps_r * laplacian;
	}

	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	}
};

// -----------------------------------------
// 9 - radiale semplice (test flusso su cerchio)
// -----------------------------------------

template <int dim>
class ExactSolution9 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		const double r = std::sqrt(x*x + y*y);

		return Ve * (r-L) / (l-L);
	};
	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	};

	virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		const double r = std::sqrt(x*x + y*y);

		Tensor<1, dim> grad;

		grad[0] = Ve / ((l-L)*r) * x;
		grad[1] = Ve / ((l-L)*r) * y;

		return grad;
};
};


template <int dim>
class RightHandSide9 : public Function<dim>{
public:
	virtual double value(const Point<dim>  &p, const unsigned int component = 0) const override{
		(void)component;
		const auto x = p[0];
		const auto y = p[1];
		const double r = std::sqrt(x*x + y*y);

		double laplacian = Ve / ((l-L)*r);
		return - eps_0 * eps_r * laplacian;
	}

	virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
		(void)component;
		AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p){
			values[p] = this->value(point_list[p]);
		}
	}
};

#endif