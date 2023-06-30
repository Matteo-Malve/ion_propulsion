#ifndef GETPOT_PRIMALSOLVER_H
#define GETPOT_PRIMALSOLVER_H

#include "Solver.h"

template <int dim>
class PrimalSolver : public Solver<dim>
{
public:
    PrimalSolver(Triangulation<dim> &       triangulation,
                 const FiniteElement<dim> & fe,
                 const Quadrature<dim> &    quadrature,
                 const Quadrature<dim - 1> &face_quadrature,
                 const Function<dim> &      rhs_function);  // tolto bdry values

    virtual void output_solution() const override;

protected:
    const SmartPointer<const Function<dim>> rhs_function;
    virtual void assemble_rhs(Vector<double> &rhs) const override;
};

// CONSTRUCTOR
template <int dim>
PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &       triangulation_,
                                const FiniteElement<dim> & fe_,
                                const Quadrature<dim> &    quadrature_,
                                const Quadrature<dim - 1> &face_quadrature_,
                                const Function<dim> &      rhs_function_)  // tolto bdry values
        : Base<dim>(triangulation_)
        , Solver<dim>(triangulation_,
                      fe_,
                      quadrature_,
                      face_quadrature_)
        , rhs_function(&rhs_function_)
{}



template <int dim>
void PrimalSolver<dim>::output_solution() const
{
    DataOut<dim> data_out;
    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "solution");
    data_out.build_patches();

    std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                      ".vtu");
    data_out.write(out, DataOutBase::vtu);
}



template <int dim>
void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const
{
    FEValues<dim> fe_values(*this->fe,
                            *this->quadrature,
                            update_values | update_quadrature_points |
                            update_JxW_values);

    const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
    const unsigned int n_q_points    = this->quadrature->size();

    Vector<double>                       cell_rhs(dofs_per_cell);
    std::vector<double>                  rhs_values(n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
        cell_rhs = 0;

        fe_values.reinit(cell);

        rhs_function->value_list(fe_values.get_quadrature_points(),
                                 rhs_values);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                                rhs_values[q_point] *               // f((x_q)
                                fe_values.JxW(q_point));            // dx

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs(local_dof_indices[i]) += cell_rhs(i);
    }
}













/*
template<int dim>
class PrimalSolver : public Solverbase<dim> {
public:
    PrimalSolver(): Solverbase<dim>(1) {};
private:

    void save_finest_mesh() override{
        cout<<"Saving Primal solver's finest mesh for later use"<<endl;
        std::ofstream out("../mesh_storage/primal_finest_mesh.vtu");
        GridOut       grid_out;
        GridOutFlags::Vtu flags(true);
        grid_out.set_flags(flags);
        if(flags.serialize_triangulation==true)
            std::cout<<" GridOutFlags::Vtu::serialize_triangulation  IS  true"<<std::endl;
        else
            std::cout<<" GridOutFlags::Vtu::serialize_triangulation  IS  false"<<std::endl;
        grid_out.write_vtu(Solverbase<dim>::triangulation, out);
        std::cout<<" Mesh written to vtu"<<endl<<endl;
        cout<<"Stored Primal Solver finest mesh for later use"<<endl;
    }
};
*/





#endif //GETPOT_PRIMALSOLVER_H
