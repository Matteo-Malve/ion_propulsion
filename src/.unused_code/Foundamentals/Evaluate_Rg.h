#ifndef ION_PROPULSION_EVALUATE_RG_H
#define ION_PROPULSION_EVALUATE_RG_H

#include "../includes&parameters_setup.h"
static GetPot redefined_6_datafile("../data_setup");

template <int dim>
class Evaluate_Rg : public Function<dim>
{
public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override{
        (void)component;

        auto x = p[0];
        auto y = p[1];
        double r = sqrt(x * x + y * y);

        double Ve = 20000;
        double Re = redefined_6_datafile("wire_radius",250e-6);

        double Rg = 0;
        if (r<2*Re)
            Rg = Ve * (2 - r/Re) * (2 - r/Re);
        return Rg;
    }
};

#endif //ION_PROPULSION_EVALUATE_RG_H
