#ifndef EVALUATION_H
#define EVALUATION_H

#include "includes.h"
#include "CSVLogger.h"


namespace IonPropulsion{
  using namespace dealii;
  namespace Evaluation{
    // ------------------------------------------------------
    // EvaluationBase
    // ------------------------------------------------------

    template <int dim>
    class EvaluationBase
    {
    public:
      virtual ~EvaluationBase() = default;

      void set_refinement_cycle(const unsigned int refinement_cycle);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                              const Vector<double> & solution) const = 0;

    protected:
      unsigned int refinement_cycle;
    };

    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------

    template <int dim>
    class PointValueEvaluation : public EvaluationBase<dim>
    {
    public:
      PointValueEvaluation(const Point<dim> &evaluation_point);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                              const Vector<double> & solution) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    private:
      const Point<dim> evaluation_point;
    };

    // ------------------------------------------------------
    // PointXDerivativeEvaluation
    // ------------------------------------------------------

    template <int dim>
    class PointXDerivativeEvaluation : public EvaluationBase<dim>
    {
    public:
      PointXDerivativeEvaluation(const Point<dim> &evaluation_point);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                              const Vector<double> & solution) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    private:
      const Point<dim> evaluation_point;
    };

    // ------------------------------------------------------
    // GridOutput
    // ------------------------------------------------------

    template <int dim>
    class GridOutput : public EvaluationBase<dim>
    {
    public:
      GridOutput(const std::string &output_name_base);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                              const Vector<double> & solution) const override;

    private:
      const std::string output_name_base;
    };


  } // namespace Evaluation
}

#endif