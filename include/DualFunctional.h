#ifndef DUALFUNCTIONAL_H
#define DUALFUNCTIONAL_H

#include "includes.h"


namespace IonPropulsion{

	using namespace dealii;

  namespace DualFunctional
  {
    // ------------------------------------------------------
    // DualFunctionalBase
    // ------------------------------------------------------
    template <int dim>
    class DualFunctionalBase : public Subscriptor
    {
    public:
      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const = 0;
    };
    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------
    template <int dim>
    class PointValueEvaluation : public DualFunctionalBase<dim>
    {
    public:
      PointValueEvaluation(const Point<dim> &evaluation_point);

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    protected:
      const Point<dim> evaluation_point;
    };
    // ------------------------------------------------------
    // PointXDerivativeEvaluation
    // ------------------------------------------------------
    template <int dim>
    class PointXDerivativeEvaluation : public DualFunctionalBase<dim>
    {
    public:
      PointXDerivativeEvaluation(const Point<dim> &evaluation_point);

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    protected:
      const Point<dim> evaluation_point;
    };

  } // namespace DualFunctional

}
#endif