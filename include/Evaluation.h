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
      EvaluationBase(const unsigned degree);
      virtual ~EvaluationBase() = default;

      void set_refinement_cycle(const unsigned int refinement_cycle);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const Vector<double> & solution,
                                                        const Triangulation<dim> &       triangulation) const = 0;

    protected:
      unsigned int refinement_cycle;
      MappingQ<dim>      mapping;

    };

    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------

    template <int dim>
    class PointValueEvaluation : public EvaluationBase<dim>
    {
    public:
      PointValueEvaluation(const unsigned degree, const Point<dim> &evaluation_point);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const Vector<double> & solution,
                                                        const Triangulation<dim> &       triangulation) const override;

      DeclException1(
        ExcEvaluationPointNotFound,
        Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    private:
      const Point<dim> evaluation_point;
    };

    // ------------------------------------------------------
    // FluxEvaluation
    // ------------------------------------------------------

    template <int dim>
    class FluxEvaluation : public EvaluationBase<dim>
    {
    public:
      FluxEvaluation(const unsigned degree, const std::set<unsigned int> &boundary_ids);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const Vector<double> & solution,
                                                        const Triangulation<dim> &       triangulation) const override;

    private:
      const std::set<unsigned int> boundary_ids;

    };


    // ------------------------------------------------------
    // GridOutput
    // ------------------------------------------------------

    template <int dim>
    class GridOutput : public EvaluationBase<dim>
    {
    public:
      GridOutput(const unsigned degree, const std::string &output_name_base);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const Vector<double> & solution,
                                                        const Triangulation<dim> &       triangulation) const override;

    private:
      const std::string output_name_base;

    };

    // ------------------------------------------------------
    // L2_error_estimate
    // ------------------------------------------------------

    template <int dim>
    class L2_error_estimate : public EvaluationBase<dim>
    {
    public:
      L2_error_estimate(const unsigned degree, const Function<dim> & analytical_solution);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const Vector<double> & solution,
                                                        const Triangulation<dim> &       triangulation) const override;

    private:
      const std::string output_name_base;
      const SmartPointer<const Function<dim>> analytical_solution;
    };

    // ------------------------------------------------------
    // H1_error_estimate
    // ------------------------------------------------------

    template <int dim>
    class H1_error_estimate : public EvaluationBase<dim>
    {
    public:
      H1_error_estimate(const unsigned degree, const Function<dim> & analytical_solution);

      virtual std::pair<std::string, double> operator()(const DoFHandler<dim> &dof_handler,
                                                        const Vector<double> & solution,
                                                        const Triangulation<dim> &       triangulation) const override;

    private:
      const std::string output_name_base;
      const SmartPointer<const Function<dim>> analytical_solution;
    };

  } // namespace Evaluation
}

#endif