#ifndef GETPOT_ERRORCONTROLLER_H
#define GETPOT_ERRORCONTROLLER_H

#include "Framework.h"
#include "PrimalSolver.h"
#include "DualSolver.h"

template <int dim>
class ErrorController : public PrimalSolver<dim>, public DualSolver<dim>
{
public:
    // CONSTRUCTOR
    ErrorController(
            Triangulation<dim> &                       coarse_grid,
        const FiniteElement<dim> &                     primal_fe,
        const FiniteElement<dim> &                     dual_fe,
        const Quadrature<dim> &                        quadrature,
        const Quadrature<dim - 1> &                    face_quadrature);
        // rhs, bdry_values, functional TO BE DEFINED inside Solvers

    // METHODS
    virtual void solve_problem() override;
    virtual unsigned int n_dofs() const override;
    virtual void refine_grid() override;
    //virtual void output_solution() const override;

private:
    // SIMPLE INTERNAL CALLS
    void solve_primal_problem();
    void solve_dual_problem();
    //
    // ... removed a lot of things
    //
    void estimate_error(Vector<float> &error_indicators) const;
    //
    // ... removed a lot of things
    //
};

template <int dim>
ErrorController<dim>::ErrorController(
        Triangulation<dim> &                           coarse_grid,
        const FiniteElement<dim> &                     primal_fe,
        const FiniteElement<dim> &                     dual_fe,
        const Quadrature<dim> &                        quadrature,
        const Quadrature<dim - 1> &                    face_quadrature,
        const Function<dim> &                          rhs_function,
        const Function<dim> &                          bv,
        const DualFunctionalBase<dim> &dual_functional)
        : Base<dim>(coarse_grid)
        , PrimalSolver<dim>(coarse_grid,
                            primal_fe,
                            quadrature,
                            face_quadrature,
                            rhs_function,
                            bv)
        , DualSolver<dim>(coarse_grid,
                          dual_fe,
                          quadrature,
                          face_quadrature,
                          dual_functional)
{}


template <int dim>
void ErrorController<dim>::solve_problem()
{
    this->solve_primal_problem();
    this->solve_dual_problem();
}


template <int dim>
void ErrorController<dim>::solve_primal_problem()
{
    PrimalSolver<dim>::solve_problem();
}

template <int dim>
void ErrorController<dim>::solve_dual_problem()
{
    DualSolver<dim>::solve_problem();
}

#endif //GETPOT_ERRORCONTROLLER_H
