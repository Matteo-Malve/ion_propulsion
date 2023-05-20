#ifndef GETPOT_DUALSOLVER_H
#define GETPOT_DUALSOLVER_H

#include "SolverBase.h"
#include "DualFunctional.h"

template<int dim>
class DualSolver : public Solverbase<dim> {
public:
    DualSolver(): Solverbase<dim>(2) {Solverbase<dim>::Nmax = 0;};
private:
    EmitterFlux<dim> J;
    void assemble_rhs() override{
        J.assemble_rhs(Solverbase<dim>::dof_handler,Solverbase<dim>::system_rhs);
    }
};




#endif //GETPOT_DUALSOLVER_H
