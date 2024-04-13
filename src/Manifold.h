#include "Functions.h"


template <int dim>
class EmitterGeometry : public ChartManifold<dim, dim, dim - 1>
{
private:
	Constants constants;
public:
	EmitterGeometry(Constants _constants): ChartManifold<dim, dim, dim - 1>(),constants(_constants){};
	virtual Point<dim - 1> pull_back(const Point<dim> &space_point) const override;

	virtual Point<dim> push_forward(const Point<dim - 1> &chart_point) const override;

	virtual std::unique_ptr<Manifold<dim, dim>> clone() const override;
};

// DEFINITIONS

template <int dim>
std::unique_ptr<Manifold<dim, dim>> EmitterGeometry<dim>::clone() const
{
	return std::make_unique<EmitterGeometry<dim>>(constants);
}

template <int dim>
Point<dim> EmitterGeometry<dim>::push_forward(const Point<dim - 1> &x) const
{
	const double y = get_emitter_height(constants.X, constants.L, constants.R, x[0]);

	Point<dim> p;
	p[0] = x[0];
	p[1] = y;

	if (dim == 3)
	{
		p[2] = x[1];
	}

	return p;
}

template <int dim>
Point<dim - 1> EmitterGeometry<dim>::pull_back(const Point<dim> &p) const
{
	Point<dim - 1> x;
	x[0] = p[0];

	if (dim == 3)
	{
		x[1] = p[2];
	}

	return x;
}