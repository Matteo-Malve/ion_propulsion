#include "includes&parameters_setup.h"

template <int dim>
class RightHandSide : public Function<dim>
{
public:
    virtual void value_list(const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {

		(void)component;

		AssertDimension (point_list.size(), values.size()); // Size check

        for (unsigned int p=0; p<point_list.size(); ++p)
        	values[p] = 0.;
	}
};


Tensor<1,2> emitter_normal(const double X, const double L, const double R, const Point<2> p) {
	//Compute the normal at point p of a circular emitter of length L
	// with circular edges of radius R centered in [-R,0] and in [-L+R,0]

	Tensor<1,2> normal;
	normal[0] = 0.;
	normal[1] = p[1];

	const double x = p[0];
	const double left_center = X - L/2. + R;
	const double right_center = X + L/2. - R;

	if (x >= right_center)
		normal[0] = x - right_center;
	else if (x <= left_center)
		normal[0] = x - left_center;

	const double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]); // The norm of this vector should always be nearly equal to R

	return normal/norm;
}



static auto evaluate_grad_Rg = [](const double X, const double L, const double R, const double Vmax,
                                    const double x, const double y) {

	Tensor<1,2> grad_Rg;

	// Gradient computed analytically by hand. Check Evaluate_Rg.h for the primitive function.
	grad_Rg[0] = 0.;
	grad_Rg[1] = 0.;

	if (x <= (X-L/2.-2.*R) || x >= (X+L/2.+2.*R) || y >= 3.*R)
		return grad_Rg;

	double d = 0.;
	const double left_center = X - L/2. + R;
	const double right_center = X + L/2. - R;

	if (x >= right_center )
		d = std::sqrt((x-right_center)*(x-right_center) + y*y);
	else if (x < right_center && x > left_center )
		d = y;
	else if (x <= left_center )
		d = std::sqrt((x-left_center)*(x-left_center) + y*y);
	else
		cout << "ERROR! Not implemented!" << endl;

	d = std::max(d,R);

	if (d < 3.*R) {
		if (x >= right_center ) {
			grad_Rg[0] = - Vmax/R * (1. - (d-R)/(2.*R)) * (x-right_center)/d; //* (1. - (d-R)/(2.*R)) * 1.5
			grad_Rg[1] = - Vmax/R * (1. - (d-R)/(2.*R)) * y/d;
		} else if (x < right_center && x > left_center ) {
			grad_Rg[0] = 0.;
			grad_Rg[1] = - Vmax/R * (1. - (d-R)/(2.*R));
		} else if (x <= left_center ) {
			grad_Rg[0] = - Vmax/R * (1. - (d-R)/(2.*R)) * (x-left_center)/d;
			grad_Rg[1] = - Vmax/R * (1. - (d-R)/(2.*R)) * y/d;
		}
	}

	  return grad_Rg;
};

double get_emitter_height(const double X, const double L, const double R, const double &p)
{
	if (p <= X - L / 2. || p >= X + L / 2.)
		return 0.;

	double y = 0;
	double x = 0;

	const double left_center = X - L / 2. + R;
	const double right_center = X + L / 2. - R;

	if (p <= left_center)
		x = p - left_center;
	else if (p >= right_center)
		x = p - right_center;

	x /= R;
	y = R * std::sqrt(1. - x * x);

	return y;
}
