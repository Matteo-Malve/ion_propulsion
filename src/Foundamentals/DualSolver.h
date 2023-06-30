#ifndef GETPOT_DUALSOLVER_H
#define GETPOT_DUALSOLVER_H

#include "DualFunctional.h"
#include "Solver.h"

template <int dim>
class DualSolver : public Solver<dim>
{
public:
    DualSolver(
            Triangulation<dim> &                           triangulation,
            const FiniteElement<dim> &                     fe,
            const Quadrature<dim> &                        quadrature,
            const Quadrature<dim - 1> &                    face_quadrature,
            const DualFunctionalBase<dim> &dual_functional);

protected:
    const SmartPointer<const DualFunctionalBase<dim>>
            dual_functional;
    virtual void assemble_rhs(Vector<double> &rhs) const override;
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
                      face_quadrature_)  // Tolto brdy condition, spostato in solver::solve()
        , dual_functional(&dual_functional_)
{}

// ASSEMBLE_rhs override
template <int dim>
void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const
{
    dual_functional->assemble_rhs(this->dof_handler, rhs);
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
