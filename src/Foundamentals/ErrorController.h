#ifndef GETPOT_ERRORCONTROLLER_H
#define GETPOT_ERRORCONTROLLER_H

#include "Framework.h"

template <int dim>
class ErrorController : public PrimalSolver<dim>, public DualSolver<dim>
{


};

#endif //GETPOT_ERRORCONTROLLER_H
