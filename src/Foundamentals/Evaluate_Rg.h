#ifndef EVALUATE_RG_H
#define EVALUATE_RG_H

#include "../includes&parameters_setup.h"


template <int dim>
class Evaluate_Rg : public Function<dim>
{
public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override{
        auto x = p[0];
        auto y = p[1];

        double r = sqrt(x * x + y * y);
        double Ve = 20000;
        double Re = 250e-6;
        double a = 1000;
        double Rg = Ve / (1 + (a*(r - Re))*(a*(r - Re)) );
        return Rg;
    }
};

#endif //EVALUATE_RG_H
