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
    using active_cell_iterator = typename DoFHandler<dim>::active_cell_iterator;
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

// CONSTRUCTOR
template <int dim>
ErrorController<dim>::ErrorController(
        Triangulation<dim> &                           coarse_grid_,
        const FiniteElement<dim> &                     primal_fe_,
        const FiniteElement<dim> &                     dual_fe_,
        const Quadrature<dim> &                        quadrature_,
        const Quadrature<dim - 1> &                    face_quadrature_,
        const Function<dim> &                          rhs_function_,        //  Not sure
        //const Function<dim> &                          bv,                  //  Not sure
        const DualFunctionalBase<dim> &dual_functional_      //  Not sure
)
        : Base<dim>(coarse_grid_)
        , PrimalSolver<dim>(coarse_grid_,
                            primal_fe_,
                            quadrature_,
                            face_quadrature_,
                            rhs_function_)
        , DualSolver<dim>(coarse_grid_,
                          dual_fe_,
                          quadrature_,
                          face_quadrature_,
                          dual_functional_)
{}

// CONSTRUCTOR CellData
template <int dim>
ErrorController<dim>::CellData::CellData(
        const SmartPointer<const FiniteElement<dim>> ptr_fe,
        const SmartPointer<const Quadrature<dim>>   ptr_quadrature,
        const SmartPointer<const Function<dim>>     ptr_right_hand_side)
        : fe_values(*ptr_fe,
                    *ptr_quadrature,
                    update_values | update_hessians | update_quadrature_points |
                    update_JxW_values)
        , right_hand_side(ptr_right_hand_side)
        , cell_residual(ptr_quadrature->size())
        , rhs_values(ptr_quadrature->size())
        , dual_weights(ptr_quadrature->size())
        , cell_laplacians(ptr_quadrature->size())
{
    cell_primal_gradients.resize(ptr_quadrature->size());
    cell_dual_gradients.resize(ptr_quadrature->size());
}

template <int dim>
ErrorController<dim>::CellData::CellData(const CellData &cell_data)
        : fe_values(cell_data.fe_values.get_fe(),
                    cell_data.fe_values.get_quadrature(),
                    update_values | update_hessians | update_quadrature_points |
                    update_JxW_values)
        , right_hand_side(cell_data.right_hand_side)
        , cell_residual(cell_data.cell_residual)
        , rhs_values(cell_data.rhs_values)
        , dual_weights(cell_data.dual_weights)
        , cell_laplacians(cell_data.cell_laplacians)
        , cell_primal_gradients(cell_data.cell_primal_gradients)
        , cell_dual_gradients(cell_data.cell_dual_gradients)

{}

template <int dim>
void ErrorController<dim>::solve_problem()
{
    this->PrimalSolver<dim>::solve_problem();
    this->DualSolver<dim>::solve_problem();
}

template <int dim>
unsigned int ErrorController<dim>::n_dofs() const
{
    return PrimalSolver<dim>::n_dofs();
}


template <int dim>
void ErrorController<dim>::refine_grid() {
    Vector<float> error_indicators(this->triangulation->n_active_cells());
    cout<<"   [ErrorController::refine_grid]Build vector to store errors"<<endl;
    estimate_error(error_indicators);
    cout<<"   [ErrorController::refine_grid]Estimate errors"<<endl;

    //compute_global_error();
    //compute_global_error_as_sum_of_cell_errors();

    for (float &error_indicator : error_indicators)
        error_indicator = std::fabs(error_indicator);

    GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                      error_indicators,
                                                      0.8,
                                                      0.02);
    this->triangulation->execute_coarsening_and_refinement();
    cout<<"   [ErrorController::refine_grid]Executed coarsening and refinement"<<endl;
}

template <int dim>
void ErrorController<dim>::estimate_error(Vector<float> &error_indicators) const
{   // INTERPOLATION
    AffineConstraints<double> dual_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler,
                                            dual_hanging_node_constraints);
    dual_hanging_node_constraints.close();
    Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs());
    cout<<"   Interpolated primal solution DOF "<<primal_solution.size()<< endl;
    FETools::interpolate(PrimalSolver<dim>::dof_handler,
                         PrimalSolver<dim>::solution,
                         DualSolver<dim>::dof_handler,
                         dual_hanging_node_constraints,
                         primal_solution);

    AffineConstraints<double> primal_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler,
                                            primal_hanging_node_constraints);
    primal_hanging_node_constraints.close();
    Vector<double> dual_weights(DualSolver<dim>::dof_handler.n_dofs());
    FETools::interpolation_difference(DualSolver<dim>::dof_handler,
                                      dual_hanging_node_constraints,
                                      DualSolver<dim>::solution,
                                      PrimalSolver<dim>::dof_handler,
                                      primal_hanging_node_constraints,
                                      dual_weights);
    cout<<"   [ErrorController::estimate_error]Interpolations done"<<endl;
    // Instance CELL DATA
    CellData cell_data(this->DualSolver<dim>::fe,
                       this->DualSolver<dim>::quadrature,
                       this->PrimalSolver<dim>::rhs_function);
    cout<<"   [ErrorController::estimate_error]Instantiated CellData"<<endl;

    // INTEGRATE OVER CELLS
    for (const auto &cell : DualSolver<dim>::dof_handler.active_cell_iterators()){
        integrate_over_cell(cell,
                            primal_solution,
                            //PrimalSolver<dim>::solution,
                            dual_weights,
                            cell_data,
                            error_indicators);
    }
    cout<<"   [ErrorController::estimate_error]Loop integral over cells done"<<endl;

}

template <int dim>
void ErrorController<dim>::integrate_over_cell(
        const active_cell_iterator &cell,
        const Vector<double> &      primal_solution,
        const Vector<double> &      dual_weights,
        CellData &                  cell_data,
        Vector<float> &             error_indicators) const
{
    cell_data.fe_values.reinit(cell);
    cell_data.right_hand_side->value_list(cell_data.fe_values.get_quadrature_points(), cell_data.rhs_values);
    cell_data.fe_values.get_function_gradients(primal_solution, cell_data.cell_primal_gradients);
    cell_data.fe_values.get_function_gradients(dual_weights, cell_data.cell_dual_gradients);


    double sum = 0;
    for (unsigned int p = 0; p < cell_data.fe_values.n_quadrature_points; ++p) {
        sum += ((cell_data.cell_primal_gradients[p]) *
                cell_data.cell_dual_gradients[p] * cell_data.fe_values.JxW(p));
    }
    error_indicators(cell->active_cell_index()) += sum;
}


#endif //GETPOT_ERRORCONTROLLER_H
