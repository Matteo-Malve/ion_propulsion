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
    protected:
      MappingQ<dim>      mapping;
    public:
      DualFunctionalBase(const unsigned mapping_degree);
      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                PETScWrappers::MPI::Vector &       rhs
                                , AffineConstraints<double> &) const = 0;

    };
    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------
    template <int dim>
    class PointValueEvaluation : public DualFunctionalBase<dim>
    {
    public:
      PointValueEvaluation(
        const unsigned mapping_degree,
        const Point<dim> &evaluation_point);

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                PETScWrappers::MPI::Vector &       rhs
                                , AffineConstraints<double> &) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    protected:
      const Point<dim> evaluation_point;
    };

    // ------------------------------------------------------
    // StandardFluxEvaluation
    // ------------------------------------------------------
    template <int dim>
    class StandardFluxEvaluation : public DualFunctionalBase<dim>
    {
    public:
      StandardFluxEvaluation(
        const unsigned mapping_degree,
        const std::set<unsigned int> &boundary_ids);

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                PETScWrappers::MPI::Vector &       rhs
                                , AffineConstraints<double> &) const override;

    protected:
      const std::set<unsigned int> boundary_ids;
    };



  } // namespace DualFunctional

}
#endif