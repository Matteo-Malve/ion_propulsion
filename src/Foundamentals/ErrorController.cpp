#include "ErrorController.h"


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

    double global_error=0.0;
    global_error=global_estimate();
    double global_error_as_sum_of_cell_errors=0.0;
    for(size_t i=0;i<error_indicators.size();i++)
        global_error_as_sum_of_cell_errors+=error_indicators[i];
    cout<<"   [ErrorController::refine_grid]Global error = "<<global_error<<endl
        <<"   [ErrorController::refine_grid]Global error as sum of cells' errors = "<<global_error_as_sum_of_cell_errors<<endl;

    for (float &error_indicator : error_indicators)
        error_indicator = std::fabs(error_indicator);
    /*
    GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                      error_indicators,
                                                      0.8,
                                                      0.02);
                                                      */
    GridRefinement::refine_and_coarsen_optimize(*this->triangulation,
                                                error_indicators,
                                                4);

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
                         primal_solution);                                      // uh projected dual FE space

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
                                      dual_weights);                            // zh-∏zh
    cout<<"   [ErrorController::estimate_error]Interpolations done"<<endl;

    // Instance CELL DATA
    CellData cell_data(this->DualSolver<dim>::fe,               // Initialized
                       this->DualSolver<dim>::quadrature,       // Empty
                       this->PrimalSolver<dim>::rhs_function);  // Empty
    cout<<"   [ErrorController::estimate_error]Instantiated CellData"<<endl;

    // INTEGRATE OVER CELLS
    for (const auto &cell : DualSolver<dim>::dof_handler.active_cell_iterators()){
        integrate_over_cell(cell,
                            primal_solution,
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
    auto & quadrature_points = cell_data.fe_values.get_quadrature_points();
    cell_data.right_hand_side->value_list(quadrature_points, cell_data.rhs_values);
    cell_data.fe_values.get_function_gradients(primal_solution, cell_data.cell_primal_gradients);
    cell_data.fe_values.get_function_gradients(dual_weights, cell_data.cell_dual_gradients);

    // Electrical permittivity of void:
    const double eps0 = 8.854*1e-12; // [F/m]

    double sum = 0;
    for (unsigned int p = 0; p < quadrature_points.size(); ++p) {
        sum +=  1.0006*eps0*1e-2*
                ((cell_data.cell_primal_gradients[p] * cell_data.cell_dual_gradients[p]  )   // Scalar product btw Tensors
                 * cell_data.fe_values.JxW(p));
    }
    error_indicators(cell->active_cell_index()) += (0 - sum);
}

template <int dim>
double ErrorController<dim>::global_estimate() const {
    // INTERPOLATION
    AffineConstraints<double> dual_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler,
                                            dual_hanging_node_constraints);
    dual_hanging_node_constraints.close();
    Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs());
    cout<<"   [ErrorController::global_estimate]Interpolated primal solution DOF "<<primal_solution.size()<< endl;
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
    cout<<"   [ErrorController::global_estimate]Interpolations done"<<endl;

    // EVALUATION
    Vector<double> temp(dual_weights.size());
    double global_error = 0.0;
    DualSolver<dim>::system_matrix.vmult(temp,dual_weights);

    if(temp.size()!=primal_solution.size())
        cout<<"PROBLEMA DIMENSIONALE"<<endl;
    for(size_t i=0;i<primal_solution.size();++i)
        global_error+=primal_solution(i)*temp(i);
    return -global_error;
}


// Instantiation

template class ErrorController<2>;



