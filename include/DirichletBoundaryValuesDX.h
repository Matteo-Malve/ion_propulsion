#ifndef ION_PROPULSION_DIRICHLETBOUNDARYVALUESDX_H
#define ION_PROPULSION_DIRICHLETBOUNDARYVALUESDX_H

#include "includes&parameters_setup.h"


template <int dim>
class DirichletBoundaryValuesDX : public Function<dim>
{
public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override{
        double h=p[1];
        return (1-(sqrt(h*h+16)/(sqrt(h*h+16)+sqrt(h*h+16*16))))*2.e4;
    }
};

#endif //ION_PROPULSION_DIRICHLETBOUNDARYVALUESDX_H
