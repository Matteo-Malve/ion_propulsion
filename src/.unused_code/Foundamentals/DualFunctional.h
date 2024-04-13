#ifndef ION_PROPULSION_DUALFUNCTIONAL_H
#define ION_PROPULSION_DUALFUNCTIONAL_H

#include "../includes&parameters_setup.h"
static GetPot redefined_3_datafile("../data_setup");

// The heart of this class. We inherit form a deal.ii class called Subscriptor
// It is a Base class containing tools for handling smart pointers
template <int dim>
class DualFunctionalBase : public Subscriptor{
public:
    virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                              Vector<double> &       rhs) const = 0;
};


// Here we build the class to use our Goal Functional as rhs.
// Evidently, we could have skipped the double inheritance, but this way we open
// ourselves to the possibility to add new Goal Functionals in the future
template <int dim>
class EmitterFlux : public DualFunctionalBase<dim>{
public:
    EmitterFlux()=default;
    virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                              Vector<double> &       rhs) const override;
};

// For details check the implementation in the corresponding .cpp file

#endif //ION_PROPULSION_DUALFUNCTIONAL_H
