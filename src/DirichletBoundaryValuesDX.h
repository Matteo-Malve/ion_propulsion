#ifndef ION_PROPULSION_DIRICHLETBOUNDARYVALUESDX_H
#define ION_PROPULSION_DIRICHLETBOUNDARYVALUESDX_H

#include "../include/includes&parameters_setup.h"


template <int dim>
class DirichletBoundaryValuesDX : public Function<dim>
{
public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
};

#endif //ION_PROPULSION_DIRICHLETBOUNDARYVALUESDX_H
