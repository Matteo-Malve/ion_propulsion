#include "PrimalSolver.h"
#include "../Evaluation.h"
#include "HelperFunctions.h"
#include "Evaluate_Rg.h"

template <int dim>
void PrimalSolver<dim>::output_solution(){

    // Only if working on main grid
    if(grid_option==1) {
        if (this->refinement_cycle == 0) {
            values.reinit(Nmax + 2);
            values(0) = 0;
        }
        values.reinit(Nmax + 2);

        // Print evaluation of V in a meaningful sample point
        Point <dim> evaluation_point(0.01, 0.002);
        Evaluation::PointValueEvaluation<dim> postprocessor(evaluation_point);
        double x_ = postprocessor(this->dof_handler, this->solution);
        std::cout << "   [PrimalSolver]Potential at sample point (" << evaluation_point[0] << "," << evaluation_point[1] << "): "
                  << std::scientific << x_ << std::defaultfloat << std::endl;
    }
    if(this->grid_option==2)
        cout<<"   No point evaluations, grid 2"<<endl;

    // WRITE SOL TO FILE.    Format: .vtu
    GradientPostprocessor<dim> gradient_postprocessor;
    DataOut <dim> data_out;
    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "Potential");
    data_out.add_data_vector(this->solution, gradient_postprocessor);
    data_out.build_patches();
    std::ofstream output("solution-" + std::to_string(this->refinement_cycle) + ".vtu");
    data_out.write_vtu(output);
}

template<int dim>
void PrimalSolver<dim>::apply_boundary_conditions() {
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    // Associate BCs to corresponding Boundary Ids. The latter are directly set in Gmsh, as physical curves
    VectorTools::interpolate_boundary_values(this->dof_handler,
                                             1, // Boundary id corresponding to Emitter
                                             Functions::ConstantFunction<dim>(0),
                                             emitter_boundary_values);

    VectorTools::interpolate_boundary_values(this->dof_handler,
                                             2,  // Boundary id corresponding to the Collector
                                             Functions::ConstantFunction<dim>(0),
                                             collector_boundary_values);
    // NOTE: We manually do the lifting.
    //       Therefore, we set the Dirichlet BCs to zero
    //       By default, deal.ii will then set all non-specified boundary points to have homogeneous Neumann BCs.
    //       That's exactly what we need.

    MatrixTools::apply_boundary_values(emitter_boundary_values,
                                       this->system_matrix,
                                       this->solution,
                                       this->system_rhs);

    MatrixTools::apply_boundary_values(collector_boundary_values,
                                       this->system_matrix,
                                       this->solution,
                                       this->system_rhs);
}

template <int dim>
void PrimalSolver<dim>::solve_problem()
{
    this->setup_system();
    this->assemble_system();                    // Assembles:  a(u0,v),  f(v),  -a(Rg,v)
    this->apply_boundary_conditions();
    this->solve_system();                       // Solves     a(u0,v) = f(v) -a(Rg,v)

    // Save uh0 for later use
    uh0 = this->solution;

    // Retrieve lifting
    Vector<double> Rg_vector(this->solution.size());
    VectorTools::interpolate(this->dof_handler,
                             Evaluate_Rg<dim>(),
                             Rg_vector);
    this->solution += Rg_vector;               // uh = u0 + Rg
}


template <int dim>
void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const {

    // Electrical permittivity of void
    const double eps0 = 8.854*1e-12; // [F/m]

    // 1) ASSEMBLE    f(v)
    FEValues <dim> fe_values(*this->fe,
                             *this->quadrature,
                             update_values | update_gradients | update_quadrature_points |
                             update_JxW_values);

    const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
    const unsigned int n_q_points = this->quadrature->size();

    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<double> rhs_values(n_q_points);
    std::vector <types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell: this->dof_handler.active_cell_iterators()) {
        cell_rhs = 0;
        fe_values.reinit(cell);
        rhs_function->value_list(fe_values.get_quadrature_points(),
                                 rhs_values);
        // Compute f_loc
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                                rhs_values[q_point] *               // f(x_q)
                                fe_values.JxW(q_point));            // dx
        // Local to Global
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs(local_dof_indices[i]) += cell_rhs(i);
        }

    // 2) ASSEMBLE   - a(Rg,v)

    // Setup
    QGauss<dim>          Rg_quadrature(2);
    FEValues<dim>        Rg_fe_values(*this->fe,
                                   Rg_quadrature,
                                   update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);
    const unsigned int Rg_n_q_points = Rg_quadrature.size();

    for (const auto &cell: this->dof_handler.active_cell_iterators()) {
        cell_rhs = 0;
        Rg_fe_values.reinit(cell);
        // Cycle over quadrature nodes
        for (unsigned int q_point = 0; q_point < Rg_n_q_points; ++q_point){
            // evaluate grad_Rg
            auto & quad_point_coords = Rg_fe_values.get_quadrature_points()[q_point];
            auto grad_Rg_xq = evaluate_grad_Rg(quad_point_coords[0],quad_point_coords[1]);
            // assemble A_loc(Rg,v)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += 1.0006*eps0*
                                (Rg_fe_values.shape_grad(i, q_point) *     // grad phi_i(x_q)
                                grad_Rg_xq *                               // grad_Rg(x_q)
                                Rg_fe_values.JxW(q_point));                // dx
        }
        // Local to Global
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs(local_dof_indices[i]) += - cell_rhs(i);
    }
}


// #######################################
// Template initialization
// #######################################
template class PrimalSolver<2>;

