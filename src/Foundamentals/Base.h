#ifndef GETPOT_BASE_H
#define GETPOT_BASE_H

#include "../includes&parameters_setup.h"
static GetPot datafile("../data_setup");

template<int dim>
class Base{
public:
    // Constructor and virtual destructor
    Base(Triangulation<dim> &coarse_grid);
    virtual ~Base() = default;

    // Pure virtual main methods
    virtual void solve_problem()                                  = 0;
    virtual void refine_grid(int step)                            = 0;
    virtual unsigned int n_dofs() const                           = 0;

    // Setter
    virtual void set_refinement_cycle(const unsigned int cycle);

    // Virtual method for output to be overwritten by PrimalSolver
    virtual void output_solution()                                = 0;

protected:
    const SmartPointer<Triangulation<dim>> triangulation;
    unsigned int refinement_cycle;

};

// CONSTRUCTOR
template <int dim>
Base<dim>::Base(Triangulation<dim> &coarse_grid_)
        : triangulation(&coarse_grid_)
        , refinement_cycle(numbers::invalid_unsigned_int)
{}

// Setter
template <int dim>
void Base<dim>::set_refinement_cycle(const unsigned int cycle)
{
    refinement_cycle = cycle;
}

#endif //GETPOT_BASE_H
