#ifndef GETPOT_PRIMALSOLVER_H
#define GETPOT_PRIMALSOLVER_H

#include "Solver.h"


template <int dim>
class PrimalSolver : public Solver<dim>
{
public:
    PrimalSolver(Triangulation<dim> &       triangulation,
                 const FiniteElement<dim> & fe,
                 const Quadrature<dim> &    quadrature,
                 const Quadrature<dim - 1> &face_quadrature,
                 const Function<dim> &      rhs_function);  // tolto bdry values
    virtual void solve_problem() override;
    virtual void output_solution() override;

protected:
    const SmartPointer<const Function<dim>> rhs_function;
    virtual void assemble_rhs(Vector<double> &rhs) const override;

private:
    Vector<float> values;
    unsigned int grid_option = datafile("Load/grid_option",1);
    unsigned int Nmax = datafile("Numerics/FEM_cycles/Nmax",10);
    double wire_radius = datafile("Mesh/wire_radius",0.025);
    const float conv_tol = datafile("Numerics/FEM_cycles/global_tolerane",1e-4);
};

// CONSTRUCTOR
template <int dim>
PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &       triangulation_,
                                const FiniteElement<dim> & fe_,
                                const Quadrature<dim> &    quadrature_,
                                const Quadrature<dim - 1> &face_quadrature_,
                                const Function<dim> &      rhs_function_)  // tolto bdry values
        : Base<dim>(triangulation_)
        , Solver<dim>(triangulation_,
                      fe_,
                      quadrature_,
                      face_quadrature_)
        , rhs_function(&rhs_function_)
{}




#endif //GETPOT_PRIMALSOLVER_H
