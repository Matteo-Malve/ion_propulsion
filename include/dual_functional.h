#ifndef DUAL_FUNCTIONAL_H
#define DUAL_FUNCTIONAL_H

#include "globals.h"
#include <deal.II/base/subscriptor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>


using namespace dealii;

template <int dim>
class DualFunctionalBase : Subscriptor {
public:
  virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                            Vector<double>        &rhs) const = 0;
};

// ------------------------------------------------------------      
// POINT evaluation
// ------------------------------------------------------------      

template <int dim>
class PointValueEvaluation : public DualFunctionalBase<dim>{
public:
  PointValueEvaluation(const Point<dim> &evaluation_point);
  virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                            Vector<double>        &rhs) const override;
  DeclException1(
    ExcEvaluationPointNotFound,
    Point<dim>,
    << "The evaluation point " << arg1
    << " was not found among the vertices of the present grid.");
protected:
  const Point<dim> evaluation_point;
};

// ------------------------------------------------------------      
// POINT dY
// ------------------------------------------------------------      

template <int dim>
class PointYDerivativeEvaluation : public DualFunctionalBase<dim>
{
public:
  PointYDerivativeEvaluation(const Point<dim> &evaluation_point);
  virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                            Vector<double>        &rhs) const override;
  DeclException1(
    ExcEvaluationPointNotFound,
    Point<dim>,
    << "The evaluation point " << arg1
    << " was not found among the vertices of the present grid.");
protected:
  const Point<dim> evaluation_point;
};

// ------------------------------------------------------------      
// AREA evaluation
// ------------------------------------------------------------  

template <int dim>
class AreaEvaluation : public DualFunctionalBase<dim>{
public:
  AreaEvaluation(const Point<dim> &center_point, const double radius);
  virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                            Vector<double>        &rhs) const override;
  DeclException1(
    ExcEvaluationPointNotFound,
    Point<dim>,
    << "The evaluation point " << arg1
    << " was not found among the vertices of the present grid.");
protected:
  const Point<dim> center_point;
  const double radius;
};

// ------------------------------------------------------------      
// BOUNDARY FLUX
// ------------------------------------------------------------      

template <int dim>
class BoundaryFluxEvaluation : public DualFunctionalBase<dim> {
public:
  BoundaryFluxEvaluation(const unsigned int boundary_id);
  virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                            Vector<double>        &rhs) const override;
protected:
  const unsigned int boundary_id;
};


#endif