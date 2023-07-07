#ifndef GETPOT_BASE_H
#define GETPOT_BASE_H

#include "../includes&parameters_setup.h"

template<int dim>
class Base{
public:
    Base(Triangulation<dim> &coarse_grid);
    virtual ~Base() = default;

    virtual void solve_problem() = 0;
    virtual void         refine_grid()                            = 0;
    virtual unsigned int n_dofs() const                           = 0;

    virtual void set_refinement_cycle(const unsigned int cycle);
    // Tolto Postprocess
    virtual void output_solution(){
        cout<<"   [Base] I don't print anything but I work"<<endl;
    }

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


template <int dim>
void Base<dim>::set_refinement_cycle(const unsigned int cycle)
{
    refinement_cycle = cycle;
}

#endif //GETPOT_BASE_H
