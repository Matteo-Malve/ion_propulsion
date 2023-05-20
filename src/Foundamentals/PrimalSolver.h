#ifndef GETPOT_PRIMALSOLVER_H
#define GETPOT_PRIMALSOLVER_H

#include "SolverBase.h"
template<int dim>
class PrimalSolver : public Solverbase<dim> {
public:
    PrimalSolver(): Solverbase<dim>(1) {};
protected:
    virtual void assemble_rhs() override{
        VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), system_matrix);
    }
    void print(){
        std::cout<<Solverbase<dim>::fe;
    }
};


#endif //GETPOT_PRIMALSOLVER_H
