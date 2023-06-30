#ifndef GETPOT_ERRORESTIMATE_H
#define GETPOT_ERRORESTIMATE_H
#include <deal.II/lac/affine_constraints.h>
#include "PrimalSolver.h"
#include "DualSolver.h"

template<int dim>
class ErrorEstimate: public PrimalSolver<dim>, DualSolver<dim>{
public:

    ErrorEstimate(
            Triangulation<dim> &                           coarse_grid,
            const FiniteElement<dim> &                     primal_fe,
            const FiniteElement<dim> &                     dual_fe,
            const Quadrature<dim> &                        quadrature,
            const Quadrature<dim - 1> &                    face_quadrature,
            const Function<dim> &                          rhs_function,
            const Function<dim> &                          boundary_values,
            const DualFunctionalBase<dim> &dual_functional):
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

    ErrorEstimate():
            uh(PrimalSolver<dim>::solution),
            z(DualSolver<dim>::solution),
            P__fe(PrimalSolver<dim>::fe),
            D__fe(DualSolver<dim>::fe),
            P__dof_handler(PrimalSolver<dim>::dof_handler),
            D__dof_handler(DualSolver<dim>::dof_handler),
            triangulation(PrimalSolver<dim>::triangulation),     // Ma devo risolvere il dual problem ad ogni passo?
            system_matrix(P__dof_handler.n_dofs(),D__dof_handler.n_dofs())
    {};

    void evaluate_err();
private:
    Vector<double> uh;
    Vector<double> z;
    double error;

    FE_Q<dim>          P__fe;
    FE_Q<dim>          D__fe;
    DoFHandler<dim>    P__dof_handler;
    DoFHandler<dim>    D__dof_handler;


    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results() const;

    Triangulation<dim> triangulation;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;


};

/*

template<int dim>
class ErrorEstimate: public PrimalSolver<dim>, DualSolver<dim>{
public:
    ErrorEstimate():
            uh(PrimalSolver<dim>::solution),
            z(DualSolver<dim>::solution),
            P__fe(PrimalSolver<dim>::fe),
            D__fe(DualSolver<dim>::fe),
            P__dof_handler(PrimalSolver<dim>::dof_handler),
            D__dof_handler(DualSolver<dim>::dof_handler),
            triangulation(PrimalSolver<dim>::triangulation),     // Ma devo risolvere il dual problem ad ogni passo?
            system_matrix(P__dof_handler.n_dofs(),D__dof_handler.n_dofs())
    {};

    void evaluate_err();
private:
    Vector<double> uh;
    Vector<double> z;
    double error;

    FE_Q<dim>          P__fe;
    FE_Q<dim>          D__fe;
    DoFHandler<dim>    P__dof_handler;
    DoFHandler<dim>    D__dof_handler;


    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results() const;

    Triangulation<dim> triangulation;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;


};

template<int dim>
auto ErrorEstimate<dim>::evaluate_err() {
    setup_system();
}

template <int dim>
void ErrorEstimate<dim>::setup_system()
{
    DynamicSparsityPattern dsp(P__dof_handler.n_dofs(),D__dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(P__dof_handler,D__dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
}

template<int dim>
void ErrorEstimate<dim>::assemble_system()
{
    QGauss<dim> quadrature_formula(P__fe.degree + 1);
    FEValues<dim> P__fe_values(P__fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values);
    FEValues<dim> D__fe_values(D__fe,
                               quadrature_formula,
                               update_values | update_gradients | update_JxW_values);

    const unsigned int P__dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int D__dofs_per_cell = fe.n_dofs_per_cell();

    FullMatrix<double> cell_matrix(P__dofs_per_cell, D__dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(P__dofs_per_cell,D__dofs_per_cell); //??

    for (const auto &cell : P__dof_handler.active_cell_iterators())
    {
        P__fe_values.reinit(cell);
        D__fe_values.reinit(cell);

        cell_matrix = 0;
        cell_rhs    = 0;

        for (const unsigned int q_index : P__fe_values.quadrature_point_indices())
        {
            for (const unsigned int i : P__fe_values.dof_indices())
                for (const unsigned int j : D__fe_values.dof_indices())
                    cell_matrix(i, j) +=
                            (P__fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                             D__fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                             fe_values.JxW(q_index));           // dx

        }
        cell->get_dof_indices(local_dof_indices);

        for (const unsigned int i : P__fe_values.dof_indices())
            for (const unsigned int j : D__fe_values.dof_indices())
                system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i, j));

        for (const unsigned int i : fe_values.dof_indices())
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<2>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}
*/

#endif //GETPOT_ERRORESTIMATE_H

