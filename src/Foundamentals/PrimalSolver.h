#ifndef GETPOT_PRIMALSOLVER_H
#define GETPOT_PRIMALSOLVER_H

#include "SolverBase.h"
template<int dim>
class PrimalSolver : public Solverbase<dim> {
public:
    PrimalSolver(): Solverbase<dim>(1) {};
private:
};


#endif //GETPOT_PRIMALSOLVER_H
