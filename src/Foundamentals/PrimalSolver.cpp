#include "PrimalSolver.h"
#include "../Evaluation.h"
#include "HelperFunctions.h"
#include "Evaluate_Rg.h"

template <int dim>
void PrimalSolver<dim>::output_solution(){
    if(grid_option==1) {
        /*
        DataOut<dim> data_out;
        data_out.attach_dof_handler(this->dof_handler);
        data_out.add_data_vector(this->solution, "solution");
        data_out.build_patches();

        std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                          ".vtu");
        data_out.write(out, DataOutBase::vtu);
        cout<<"Output solution done\n";
        */

        //-------------------------
        if (this->refinement_cycle == 0) {
            values.reinit(Nmax + 2);
            values(0) = 0;
        }
        values.reinit(Nmax + 2);

        Point <dim> evaluation_point(0.01, 0.002);

        // Print sample V
        Evaluation::PointValueEvaluation<dim> postprocessor(evaluation_point);
        double x_ = postprocessor(this->dof_handler, this->solution);
        std::cout << "   [PrimalSolver]Potential at (" << evaluation_point[0] << "," << evaluation_point[1] << "): "
                  << std::scientific << x_ << std::defaultfloat << std::endl;

        // Upon reaching of convergence, print results
        if (this->refinement_cycle == Nmax) {

            Point <dim> sample(wire_radius, 0.);
            Tensor<1, dim> E = VectorTools::point_gradient(this->dof_handler, this->solution, sample);
            std::cout << "   [PrimalSolver]Electric field in (" << sample[0] << "," << sample[1] << "): " << -E
                      << ", magnitude: " << L2Norm(E) << std::endl;
        }
    }
    if(this->grid_option==2)
        cout<<"   No point evaluations, grid 2"<<endl;

    // WRITE SOL TO .vtu
    GradientPostprocessor<dim> gradient_postprocessor;
    DataOut <dim> data_out;
    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "Potential");
    data_out.add_data_vector(this->solution, gradient_postprocessor);
    data_out.build_patches();
    std::ofstream output("solution-" + std::to_string(this->refinement_cycle) + ".vtu");
    data_out.write_vtu(output);
}



template <int dim>
void PrimalSolver<dim>::solve_problem()
{
    this->setup_system();
    this->assemble_system(); // Assemble a(u0,v) f(v) -a(Rg,v)

    this->apply_boundary_conditions(); // Applica a u0 BC Neumann: sotto, sx e sopra
    this->solve_system();   // Solve     a(u0,v) = f(v) -a(Rg,v)

    // Retrieve lifting
    Vector<double> Rg_vector(this->solution.size());
    VectorTools::interpolate(this->dof_handler,
                             Evaluate_Rg<dim>(),
                             Rg_vector);
    this->solution += Rg_vector;                 // uh = u0 + Rg     Si perdono le BC di Neumann essendo Rg radiale
}
/*
auto evaluate_Rg = [](double x, double y) {
    double r = sqrt(x * x + y * y);
    double Ve = 20000;
    double Re = 250e-6;
    double a = 100;
    double Rg = Ve / (1 + (a*(r - Re))*(a*(r - Re)) );
    return Rg;
};*/

auto evaluate_grad_Rg = [](double x, double y) {
    double r = sqrt(x * x + y * y);
    Tensor<1,2> grad_Rg;
    double Ve = 20000;
    double Re = 250e-6;
    double a = 100000;
    double dfdr = - Ve / ( (1 +  (a*(r - Re))*(a*(r - Re)) )*(1+ (a*(r - Re))*(a*(r - Re)) ) ) *  a*a*2*(r-Re);
    grad_Rg[0] = dfdr * x / r;
    grad_Rg[1] = dfdr * y / r;
    return grad_Rg;
};

template <int dim>
void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const {
    const double eps0 = 8.854*1e-12; // [F/m]
    // f(v)

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


        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) +=(fe_values.shape_value(i, q_point) * // phi_i(x_q)
                                rhs_values[q_point] *               // f((x_q)
                                fe_values.JxW(q_point));            // dx

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs(local_dof_indices[i]) += cell_rhs(i);
    }

    // -----------------------------

    // -a(Rg,v)

    QGauss<dim>          Rg_quadrature(4);
    FEValues<dim>        Rg_fe_values(*this->fe,
                                   Rg_quadrature,
                                   update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);

    const unsigned int Rg_n_q_points = Rg_quadrature.size();

    for (const auto &cell: this->dof_handler.active_cell_iterators()) {
        cell_rhs = 0;
        Rg_fe_values.reinit(cell);

        // A_loc(Rg,v)
        for (unsigned int q_point = 0; q_point < Rg_n_q_points; ++q_point){
            // evaluate evaluate_grad_Rg
            auto & quad_point_coords = Rg_fe_values.get_quadrature_points()[q_point];
            auto grad_Rg_xq = evaluate_grad_Rg(quad_point_coords[0],quad_point_coords[1]);
            // assemble a(Rg,v)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += 1.0006*eps0*
                                (Rg_fe_values.shape_grad(i, q_point) *      // grad phi_i(x_q)
                                grad_Rg_xq *                               // grad_Rg(x_q)
                                Rg_fe_values.JxW(q_point));                // dx
        }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs(local_dof_indices[i]) += - cell_rhs(i);
    }

    // Prepare Rg_vector for later use
/*
    Vector<double> cell_Rg_vector(9);

    for (const auto &cell: this->dof_handler.active_cell_iterators()) {
        cell_Rg_vector = 0;



            auto & quad_point_coords = Rg_fe_values.get_quadrature_points()[q_point];
            auto cell_Rg_vector[q_point] = evaluate_Rg(quad_point_coords[0],quad_point_coords[1]);

        cell->get_dof_indices(local_dof_index);
            this->Rg_vector(local_dof_index) =



    }

    for (unsigned int q_point = 0; q_point < Rg_n_q_points; ++q_point){
        // evaluate evaluate_grad_Rg
        auto & quad_point_coords = fe_values.get_quadrature_points()[q_point];
        auto Rg_xq = evaluate_Rg(quad_point_coords[0],quad_point_coords[1]);
        // assemble a(Rg,v)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) += (Rg_fe_values.shape_value(i, q_point) *      // grad phi_i(x_q)
                            Rg_xq *                               // grad_Rg(x_q)
                            Rg_fe_values.JxW(q_point));                // dx
    }

    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
        rhs(local_dof_indices[i]) += - cell_rhs(i);
}
*/
}


// #######################################
// Template initialization
// #######################################
template class PrimalSolver<2>;

