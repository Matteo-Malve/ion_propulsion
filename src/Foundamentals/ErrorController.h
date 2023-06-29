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
        const Quadrature<dim - 1> &                    face_quadrature,
        const Function<dim> &                          rhs_function,
        const DualFunctionalBase<dim> &dual_functional);
    // Bdry values to be done in solvers

    // METHODS
    virtual void solve_problem() override;
    virtual unsigned int n_dofs() const override;
    virtual void refine_grid() override;
    //virtual void output_solution() const override;                // TO BE DONE

private:
    // SIMPLE INTERNAL CALLS
    void solve_primal_problem();
    void solve_dual_problem();
    //
    // ... removed a lot of things
    //
    void estimate_error(Vector<float> &error_indicators) const;     // TO BE DONE
    //
    // ... removed a lot of things
    //
};

// CONSTRUCTOR
template <int dim>
ErrorController<dim>::ErrorController(
        Triangulation<dim> &                           coarse_grid,
        const FiniteElement<dim> &                     primal_fe,
        const FiniteElement<dim> &                     dual_fe,
        const Quadrature<dim> &                        quadrature,
        const Quadrature<dim - 1> &                    face_quadrature,
        const Function<dim> &                          rhs_function,        //  Not sure
        //const Function<dim> &                          bv,                  //  Not sure
        const DualFunctionalBase<dim> &dual_functional      //  Not sure
        )
        : Base<dim>(coarse_grid)
        , PrimalSolver<dim>(coarse_grid,
                            primal_fe,
                            quadrature,
                            face_quadrature,
                            rhs_function)
        , DualSolver<dim>(coarse_grid,
                          dual_fe,
                          quadrature,
                          face_quadrature,
                          dual_functional)
{}

template <int dim>
void ErrorController<dim>::solve_problem()
{
    this->PrimalSolver<dim>::solve_problem();
    this->DualSolver<dim>::solve_problem();
}

/*
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
*/
template <int dim>
unsigned int ErrorController<dim>::n_dofs() const
{
    return PrimalSolver<dim>::n_dofs();
}


template <int dim>
void ErrorController<dim>::refine_grid() {
    Vector<float> error_indicators(this->triangulation->n_active_cells());
    estimate_error(error_indicators);

    for (float &error_indicator : error_indicators)
        error_indicator = std::fabs(error_indicator);

    GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                      error_indicators,
                                                      0.8,
                                                      0.02);
    this->triangulation->execute_coarsening_and_refinement();
}

template <int dim>
void ErrorController<dim>::estimate_error(Vector<float> &error_indicators) const{
    // TO BE DONE
}

#endif //GETPOT_ERRORCONTROLLER_H
