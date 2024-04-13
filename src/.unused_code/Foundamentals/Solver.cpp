#include "Solver.h"

template <int dim>
Solver<dim>::~Solver()
{
    dof_handler.clear();
}

template <int dim>
unsigned int Solver<dim>::n_dofs() const
{
    return dof_handler.n_dofs();
}

template <int dim>
void Solver<dim>::setup_system()
{
    dof_handler.distribute_dofs(*fe);
    std::cout << "   [Solver]Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void Solver<dim>::assemble_system()
{
    FEValues<dim> fe_values(*this->fe,
                            *this->quadrature,
                            update_values | update_gradients | update_JxW_values);
    const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        fe_values.reinit(cell);
        cell_matrix = 0;

        // Electrical permittivity of void:
        const double eps0 = 8.854*1e-12; // [F/m]

        // Compute A_loc
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
            for (const unsigned int i : fe_values.dof_indices())
                for (const unsigned int j : fe_values.dof_indices())
                    cell_matrix(i, j) +=
                            1.0006*eps0*
                            (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                             fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                             fe_values.JxW(q_index));           // dx
        }
        cell->get_dof_indices(local_dof_indices);

        // Add contribution to A_global
        for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
                system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i, j));
    }

    // Call assemble_rhs (overwritten by Primal and Dual solvers)
    this->assemble_rhs(this->system_rhs);

    cout<<"   [Solver::assemble_system]Assembled the rhs"<<endl;
}


template <int dim>
void Solver<dim>::solve_system()
{
    // Build CG Solver
    SolverControl            solver_control(2000, 1e-6 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(solver_control);

    // Construct SSOR preconditioner
    double relaxation_parameter = 1.5;
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, relaxation_parameter);

    // Solve Linear System with CG and SSOR preconditioner
    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    // Print iterations
    cout<<"   [Solver]"<<solver_control.last_step()
        <<" CG iterations needed to obtain convergence."<<endl;
}




// #######################################
// Template initialization
// #######################################
template class Solver<2>;

