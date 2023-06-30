#include "Solver.h"



template<int dim>
void Solver<dim>::apply_boundary_conditions() {
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    VectorTools::interpolate_boundary_values(dof_handler,
                                             1, // Boundary corrispondente all'emettitore, definito sopra
                                             Functions::ConstantFunction<dim>(2.e+4), // Valore di potenziale all'emettitore (20 kV)
                                             emitter_boundary_values);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             2,  // Boundary corrispondente al collettore, definito sopra
                                             Functions::ConstantFunction<dim>(1.6e+4), // Valore di potenziale al collettore (0 V)
            //DirichletBoundaryValuesDX<dim>(),
                                             collector_boundary_values);

    /* Le condizioni sopra sono condizioni di Dirichlet
    Su il restante bordo del dominio, dove non è imposta esplicitamente nessuna condizione,
    viene imposta automaticamente da deal.II una condizione di Neumann: gradiente nullo */

    MatrixTools::apply_boundary_values(emitter_boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);

    MatrixTools::apply_boundary_values(collector_boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}


template <int dim>
Solver<dim>::~Solver()
{
    dof_handler.clear();
}


template <int dim>
void Solver<dim>::solve_problem()
{
    this->apply_boundary_conditions();      // OURS

    setup_system();
    assemble_system();
    solve_system();
    //output_results();
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
    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
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
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        fe_values.reinit(cell);

        cell_matrix = 0;
        cell_rhs    = 0;

        for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
            for (const unsigned int i : fe_values.dof_indices())
                for (const unsigned int j : fe_values.dof_indices())
                    cell_matrix(i, j) +=
                            (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                             fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                             fe_values.JxW(q_index));           // dx

            for (const unsigned int i : fe_values.dof_indices())
                cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                                1. *                                // f(x_q)
                                fe_values.JxW(q_index));            // dx
        }
        cell->get_dof_indices(local_dof_indices);

        for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
                system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i, j));

        for (const unsigned int i : fe_values.dof_indices())
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}


template <int dim>
void Solver<dim>::solve_system()
{
    SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}




// #######################################
// Template initialization
// #######################################
template class Solver<2>;

