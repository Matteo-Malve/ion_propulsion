#include "DirichletBoundaryValuesDX.h"

template <int dim>
double DirichletBoundaryValuesDX<dim>::value(const Point<dim> &p) const{
    double h=p[1];
    return (1-(sqrt(h*h+16)/(sqrt(h*h+16)+sqrt(h*h+16*16))))*2.e4;
}

// #######################################
// Template initialization
// #######################################
template RightDirichletBoundaryValues<2>;