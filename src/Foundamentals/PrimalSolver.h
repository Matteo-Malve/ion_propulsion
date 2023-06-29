#ifndef GETPOT_PRIMALSOLVER_H
#define GETPOT_PRIMALSOLVER_H

#include "SolverBase.h"
#include "Base.h"

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






#endif //GETPOT_PRIMALSOLVER_H
