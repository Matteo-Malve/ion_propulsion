#ifndef ION_PROPULSION_PRIMALSOLVER_H
#define ION_PROPULSION_PRIMALSOLVER_H

#include "Solver.h"


template <int dim>
class PrimalSolver : public Solver<dim>
{
public:
    PrimalSolver(Triangulation<dim> &       triangulation,
                 const FiniteElement<dim> & fe,
                 const Quadrature<dim> &    quadrature,
                 const Quadrature<dim - 1> &face_quadrature,
                 const Function<dim> &      rhs_function);
    virtual void solve_problem() override;
    virtual void output_solution() override;

protected:
    const SmartPointer<const Function<dim>> rhs_function;
    virtual void assemble_rhs(Vector<double> &rhs) const override;
    virtual void apply_boundary_conditions() override;
    Vector<double>   uh0;
private:
    Vector<float> values;
    unsigned int grid_option = datafile("grid_option",1);
    unsigned int Nmax = datafile("Nmax",10);


};

// CONSTRUCTOR
template <int dim>
PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &       triangulation_,
                                const FiniteElement<dim> & fe_,
                                const Quadrature<dim> &    quadrature_,
                                const Quadrature<dim - 1> &face_quadrature_,
                                const Function<dim> &      rhs_function_)
        : Base<dim>(triangulation_)
        , Solver<dim>(triangulation_,
                      fe_,
                      quadrature_,
                      face_quadrature_)
        , rhs_function(&rhs_function_)
{}




#endif //ION_PROPULSION_PRIMALSOLVER_H
