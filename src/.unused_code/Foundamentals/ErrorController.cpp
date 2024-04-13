#include "ErrorController.h"
#include "Evaluate_Rg.h"

// CONSTRUCTOR
template <int dim>
ErrorController<dim>::ErrorController(
        Triangulation<dim> &                           coarse_grid_,
        const FiniteElement<dim> &                     primal_fe_,
        const FiniteElement<dim> &                     dual_fe_,
        const Quadrature<dim> &                        quadrature_,
        const Quadrature<dim - 1> &                    face_quadrature_,
        const Function<dim> &                          rhs_function_,
        const DualFunctionalBase<dim> &dual_functional_
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

// COPY CONSTRUCTOR CellData
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

// Methods definitions:

template <int dim>
void ErrorController<dim>::solve_problem()
{
    cout<<"   [ErrorController::solve_problem] Begin solving Primal Problem"<<endl;
    this->PrimalSolver<dim>::solve_problem();
    cout<<"   [ErrorController::solve_problem] Begin solving Dual Problem"<<endl;
    this->DualSolver<dim>::solve_problem();
}

template <int dim>
unsigned int ErrorController<dim>::n_dofs() const
{
    return PrimalSolver<dim>::n_dofs();
}


template <int dim>
void ErrorController<dim>::refine_grid(unsigned int algorithm) {

    // Prepare vector to store error values
    Vector<float> error_indicators(this->triangulation->n_active_cells());
    // Compute local errors
    estimate_error(error_indicators);

    // Compute global error estimate
    double global_error=0.0;
    global_error=global_estimate();

    // Sum contribution of each cell's local error to get a global estimate
    double global_error_as_sum_of_cell_errors=0.0;
    for(size_t i=0;i<error_indicators.size();i++)
        global_error_as_sum_of_cell_errors+=error_indicators[i];
    global_error_as_sum_of_cell_errors=abs(global_error_as_sum_of_cell_errors);

    // Output the two derived global estimates
    cout<<"   [ErrorController::refine_grid]Global error = "<<global_error<<endl
        <<"   [ErrorController::refine_grid]Global error as sum of cells' errors = "<<global_error_as_sum_of_cell_errors<<endl;

    // Take absolute value of each error
    for (float &error_indicator : error_indicators)
        error_indicator = std::fabs(error_indicator);

    if(algorithm>1.5)
        GridRefinement::refine_and_coarsen_optimize(*this->triangulation,
                                                    error_indicators,
                                                    4);
    else
        GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                        error_indicators,
                                                        0.8,
                                                        0.02);

    // Execute refinement
    this->triangulation->execute_coarsening_and_refinement();
    cout<<"   [ErrorController::refine_grid]Executed coarsening and refinement"<<endl;
}



template <int dim>
void ErrorController<dim>::estimate_error(Vector<float> &error_indicators) const
{
    // INTERPOLATION

    // Project uh0 on dual solution FE space
    AffineConstraints<double> dual_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler,
                                            dual_hanging_node_constraints);
    dual_hanging_node_constraints.close();
    Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs());

    FETools::interpolate(PrimalSolver<dim>::dof_handler,
                         PrimalSolver<dim>::uh0,
                         DualSolver<dim>::dof_handler,
                         dual_hanging_node_constraints,
                         primal_solution);

    // Subtract from dual solution its projection on the primal solution FE space
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

    // Instance CELL DATA
    CellData cell_data(this->DualSolver<dim>::fe,               // Initialized
                       this->DualSolver<dim>::quadrature,       // Empty
                       this->PrimalSolver<dim>::rhs_function);  // Empty

    // Compute Rg_plus_uh0hat
    Vector<double> Rg_plus_uh0hat(primal_solution.size());
    Vector<double> Rg_vector(primal_solution.size());
    VectorTools::interpolate(DualSolver<dim>::dof_handler,
                             Evaluate_Rg<dim>(),
                             Rg_vector);
    for(size_t i=0;i<primal_solution.size();i++)
        Rg_plus_uh0hat(i) = Rg_vector(i) + primal_solution(i);

    // INTEGRATE OVER CELLS
    for (const auto &cell : DualSolver<dim>::dof_handler.active_cell_iterators()){
        integrate_over_cell(cell,
                            Rg_plus_uh0hat,
                            dual_weights,
                            cell_data,
                            error_indicators);
    }

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

    // Numerically approximate the integral of the scalar product between the gradients of the two
    double sum = 0;
    for (unsigned int p = 0; p < quadrature_points.size(); ++p) {
        sum +=  1.0006*eps0*
                ((cell_data.cell_primal_gradients[p] * cell_data.cell_dual_gradients[p]  )   // Scalar product btw Tensors
                 * cell_data.fe_values.JxW(p));
    }
    error_indicators(cell->active_cell_index()) += (0 - sum);
}

template <int dim>
double ErrorController<dim>::global_estimate() const {

    // INTERPOLATION
    // Project uh0 on dual solution FE space
    AffineConstraints<double> dual_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler,
                                            dual_hanging_node_constraints);
    dual_hanging_node_constraints.close();
    Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs());
    cout<<"   [ErrorController::global_estimate]Interpolated primal solution DOF "<<primal_solution.size()<< endl;
    FETools::interpolate(PrimalSolver<dim>::dof_handler,
                         PrimalSolver<dim>::uh0,
                         DualSolver<dim>::dof_handler,
                         dual_hanging_node_constraints,
                         primal_solution);

    // Subtract from dual solution its projection on the primal solution FE space
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

    // EVALUATION
    Vector<double> temp(dual_weights.size());
    double dx = 0.0; double lx = 0.0;

    // Compute   u' A z                             NOTE: z is actually the dual weights (zh-∏zh)
    DualSolver<dim>::system_matrix.vmult(temp,dual_weights);    // A z
    if(temp.size()!=primal_solution.size())
        cout<<"DIMENSIONALITY PROBLEM"<<endl;
    for(size_t i=0;i<primal_solution.size();++i)
        dx+=primal_solution(i)*temp(i);           // u' A z

    Vector<double> F(dual_weights.size());

    // Compute   a(Rg,φ)

    // NOTE: Almost identical to computation of a(Rg,v) in Primal Solver's rhs
    //       Only here we use shape functions from the dual FE space
    const double eps0 = 8.854*1e-12; // [F/m]
    const unsigned int dofs_per_cell = DualSolver<dim>::dof_handler.get_fe().n_dofs_per_cell();
    Vector<double> cell_F(dofs_per_cell);
    std::vector <types::global_dof_index> local_dof_indices(dofs_per_cell);
    QGauss<dim>          Rg_quadrature(5);
    FEValues<dim>        Rg_fe_values(DualSolver<dim>::dof_handler.get_fe(),
                                      Rg_quadrature,
                                      update_values | update_gradients | update_quadrature_points |
                                      update_JxW_values);
    const unsigned int Rg_n_q_points = Rg_quadrature.size();
    for (const auto &cell: DualSolver<dim>::dof_handler.active_cell_iterators()) {
        cell_F = 0;
        Rg_fe_values.reinit(cell);
        for (unsigned int q_point = 0; q_point < Rg_n_q_points; ++q_point){
            // evaluate_grad_Rg
            auto & quad_point_coords = Rg_fe_values.get_quadrature_points()[q_point];
            auto grad_Rg_xq = evaluate_grad_Rg(quad_point_coords[0],quad_point_coords[1]);
            // assemble a(Rg,φ)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_F(i) += 1.0006*eps0*
                               (Rg_fe_values.shape_grad(i, q_point) *      // grad phi_i(x_q)
                                grad_Rg_xq *                               // grad_Rg(x_q)
                                Rg_fe_values.JxW(q_point));                // dx
        }
        // Local to Global
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            F(local_dof_indices[i]) += cell_F(i);
    }

    // Compute     z' F    or    z•F
    for(size_t i=0;i<dual_weights.size();i++)
        lx += dual_weights(i) * F(i);


    // Return             η = |r(z)|  = | - z' F - u' A z |
    // or, precisely,     η = |r(zk-∏zk)|  = | - (zk-∏zk)' F - uh0' A (zk-∏zk) |
    return abs(-lx -dx);
}




// #######################################
// Template initialization
// #######################################
template class ErrorController<2>;



