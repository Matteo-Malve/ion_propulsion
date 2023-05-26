#ifndef GETPOT_DUALSOLVER_H
#define GETPOT_DUALSOLVER_H

#include "SolverBase.h"
#include "DualFunctional.h"

template<int dim>
class DualSolver : public Solverbase<dim> {
public:
    DualSolver(): Solverbase<dim>(2) {Solverbase<dim>::Nmax = 0;};
    void run() override;
private:
    EmitterFlux<dim> J;
    void assemble_rhs() override{
        J.assemble_rhs(Solverbase<dim>::dof_handler,Solverbase<dim>::system_rhs);
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
    Solverbase<dim>::system_matrix.copy_from(Solverbase<dim>::laplace_matrix);
    Solverbase<dim>::system_matrix *= 1.0006*eps0*1e-3;
    Solverbase<dim>::constraints.condense(Solverbase<dim>::system_matrix, Solverbase<dim>::system_rhs);
    Solverbase<dim>::apply_boundary_conditions();


    std::cout << std::endl
            << "   Number of active cells:       " << Solverbase<dim>::finest_mesh.n_active_cells() << std::endl
            << "   Number of degrees of freedom: " << Solverbase<dim>::dof_handler.n_dofs() << std::endl;

    Solverbase<dim>::solve();
    Solverbase<dim>::output_results();

    std::cout << "   Elapsed CPU time: " << Solverbase<dim>::timer.cpu_time() << " seconds.\n";

}



#endif //GETPOT_DUALSOLVER_H
