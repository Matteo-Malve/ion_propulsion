#ifndef ION_PROPULSION_ERRORCONTROLLER_H
#define ION_PROPULSION_ERRORCONTROLLER_H

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

    using active_cell_iterator = typename DoFHandler<dim>::active_cell_iterator;
    // METHODS
    virtual void solve_problem() override;
    virtual unsigned int n_dofs() const override;
    virtual void refine_grid(unsigned int algorithm) override;
    virtual void output_solution() override;
private:
    void estimate_error(Vector<float> &error_indicators) const;
    double global_estimate() const;
    struct CellData
    {
        FEValues<dim>                           fe_values;
        const SmartPointer<const Function<dim>> right_hand_side;

        std::vector<double> cell_residual;
        std::vector<double> rhs_values;
        std::vector<double> dual_weights;
        std::vector<double> cell_laplacians;    //useless
        typename std::vector<Tensor<1, dim>> cell_primal_gradients;
        typename std::vector<Tensor<1, dim>> cell_dual_gradients;
        CellData(const SmartPointer<const FiniteElement<dim>> ptr_fe,
                 const SmartPointer<const Quadrature<dim>>   ptr_quadrature,
                 const SmartPointer<const Function<dim>>     ptr_right_hand_side);
        CellData(const CellData &cell_data);
    };
    void integrate_over_cell(
            const active_cell_iterator &cell,
            const Vector<double> &      primal_solution,
            const Vector<double> &      dual_weights,
            CellData &                  cell_data,
            Vector<float> &             error_indicators) const;
};

template <int dim>
void ErrorController<dim>::output_solution()
{
    this->PrimalSolver<dim>::output_solution();
    this->DualSolver<dim>::output_solution();
}

#endif //ION_PROPULSION_ERRORCONTROLLER_H
