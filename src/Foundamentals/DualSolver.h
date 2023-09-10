#ifndef GETPOT_DUALSOLVER_H
#define GETPOT_DUALSOLVER_H

#include "DualFunctional.h"
#include "Solver.h"
#include "HelperFunctions.h"

template <int dim>
class DualSolver : public Solver<dim>
{
public:
    // Constructor
    DualSolver(
            Triangulation<dim> &                           triangulation,
            const FiniteElement<dim> &                     fe,
            const Quadrature<dim> &                        quadrature,
            const Quadrature<dim - 1> &                    face_quadrature,
            const DualFunctionalBase<dim> &dual_functional);
    // Override solve method
    virtual void solve_problem() override;
    virtual void output_solution() override;
protected:
    // DualFunctional class comes into play
    const SmartPointer<const DualFunctionalBase<dim>>
            dual_functional;
    virtual void assemble_rhs(Vector<double> &rhs) const override;
    virtual void apply_boundary_conditions() override;
};

// CONSTRUCTOR
template <int dim>
DualSolver<dim>::DualSolver(
        Triangulation<dim> &                           triangulation_,
        const FiniteElement<dim> &                     fe_,
        const Quadrature<dim> &                        quadrature_,
        const Quadrature<dim - 1> &                    face_quadrature_,
        const DualFunctionalBase<dim> &dual_functional_)
        : Base<dim>(triangulation_)
        , Solver<dim>(triangulation_,
                      fe_,
                      quadrature_,
                      face_quadrature_)
        , dual_functional(&dual_functional_)
{}

// ASSEMBLE_rhs override
template <int dim>
void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const
{
    dual_functional->assemble_rhs(this->dof_handler, rhs);

}

template<int dim>
void DualSolver<dim>::apply_boundary_conditions() {
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    VectorTools::interpolate_boundary_values(this->dof_handler,
                                             1, // Boundary corrispondente all'emettitore, definito sopra
                                             Functions::ConstantFunction<dim>(0), // Valore di potenziale all'emettitore (20 kV)
                                             emitter_boundary_values);

    VectorTools::interpolate_boundary_values(this->dof_handler,
                                             2,  // Boundary corrispondente al collettore, definito sopra
                                             Functions::ConstantFunction<dim>(0), // Valore di potenziale al collettore (0 V)
            //DirichletBoundaryValuesDX<dim>(),
                                             collector_boundary_values);

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
void DualSolver<dim>::solve_problem()
{
    this->setup_system();
    this->assemble_system();


    int flag = 0;
    for(int i=0;i<this->system_rhs.size();i++)
        if(this->system_rhs(i)!=0)
            flag++;
    cout<<"Nonzero elements of rhs BEFORE apply_boundary_conditions() = "<<flag<<endl;


    this->apply_boundary_conditions();

    flag = 0;
    for(int i=0;i<this->system_rhs.size();i++)
        if(this->system_rhs(i)!=0)
            flag++;
    cout<<"Nonzero elements of rhs AFTER apply_boundary_conditions() = "<<flag<<endl;



    //for(int i=0;i<100;i++)
    //    cout<<this->system_rhs(i)<<endl;

    this->solve_system();

}

template <int dim>
void DualSolver<dim>::output_solution()
{

    // WRITE SOL TO .vtu
    GradientPostprocessor<dim> gradient_postprocessor;
    DataOut <dim> data_out;
    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "Potential");
    data_out.add_data_vector(this->solution, gradient_postprocessor);
    data_out.build_patches();
    std::ofstream output("dual_solution-" + std::to_string(this->refinement_cycle) + ".vtu");
    data_out.write_vtu(output);

    std::ofstream out("../mesh_storage/dual_mesh_cycle_" + std::to_string(this->refinement_cycle) + ".vtu");
    GridOut       grid_out;
    GridOutFlags::Vtu flags(true);
    grid_out.set_flags(flags);
    grid_out.write_vtu(*this->triangulation, out);
    std::cout<<" Mesh written to vtu"<<endl<<endl;

}








/*
template<int dim>
class DualSolver : public Solverbase<dim> {
public:
    DualSolver(): Solverbase<dim>(2) {Solverbase<dim>::Nmax = 0;};
    void run() override;

    using Solverbase<dim>::dof_handler;
    using Solverbase<dim>::system_rhs;
    using Solverbase<dim>::laplace_matrix;
    using Solverbase<dim>::system_matrix;
    using Solverbase<dim>::finest_mesh;
    using Solverbase<dim>::timer;
    using Solverbase<dim>::constraints;
    using Solverbase<dim>::apply_boundary_conditions;
    using Solverbase<dim>::solve;
    using Solverbase<dim>::output_results;

private:
    EmitterFlux<dim> J;

    void assemble_rhs() override{
        J.assemble_rhs(dof_handler,system_rhs);
    }
};

template<int dim>
void DualSolver<dim>::run()
{
    // Load finest Mesh from Primal Solver
    cout<<"Looking for the stored Finest Mesh from Primal Solver"<<endl;
    struct stat sb;
    int is_present = stat("../mesh_storage/primal_finest_mesh.vtu",&sb);
    if(is_present==-1) {
        std::cerr << " Mesh NOT found!" << std::endl;
    }
    else if(is_present==0){
        std::cerr<<" Mesh found"<<std::endl<<" Prepare import"<<endl;
        std::ifstream input_file("../mesh_storage/primal_finest_mesh.vtu");
        GridIn<dim>       grid_in;
        grid_in.attach_triangulation(Solverbase<dim>::triangulation);
        grid_in.read_vtu(input_file);
        std::cerr<<" Grid imported"<<std::endl;
    } else
        std::cerr<<" File not found nor not found. Anomaly."<<endl;


    Solverbase<dim>::setup_system();
    assemble_rhs();
    const double eps0 = 8.854*1e-12; // [F/m]
    system_matrix.copy_from(laplace_matrix);
    system_matrix *= 1.0006*eps0*1e-3;
    constraints.condense(system_matrix, system_rhs);
    apply_boundary_conditions();


    std::cout << std::endl
            << "   Number of active cells:       " << finest_mesh.n_active_cells() << std::endl
            << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

    solve();
    output_results();

    std::cout << "   Elapsed CPU time: " << timer.cpu_time() << " seconds.\n";

}

*/

#endif //GETPOT_DUALSOLVER_H
