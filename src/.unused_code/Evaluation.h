#ifndef ION_PROPULSION_EVALUATION_H
#define ION_PROPULSION_EVALUATION_H

#include "includes&parameters_setup.h"

using namespace dealii;

namespace Evaluation
{
    template <int dim>
    class EvaluationBase
    {
    public:
        virtual ~EvaluationBase() = default;
        virtual void operator()(const DoFHandler<dim> &dof_handler,
                                const Vector<double> & solution) const = 0;


    protected:
        unsigned int refinement_cycle;
    };

    template <int dim>
    class PointValueEvaluation : public EvaluationBase<dim>
    {
    public:
        PointValueEvaluation(const Point<dim> &evaluation_point);

        virtual void operator()(const DoFHandler<dim> &dof_handler,
                                const Vector<double> & solution) const override;
        double operator()(const DoFHandler<dim> &dof_handler,
                          const Vector<double> & solution);
        DeclException1(
                ExcEvaluationPointNotFound,
                Point<dim>,
        << "The evaluation point " << arg1
        << " was not found among the vertices of the present grid.");

    private:
        const Point<dim> evaluation_point;
    };


    template <int dim>
    PointValueEvaluation<dim>::PointValueEvaluation(
            const Point<dim> &evaluation_point)
            : evaluation_point(evaluation_point)
    {}



    template <int dim>
    void
    PointValueEvaluation<dim>::operator()(const DoFHandler<dim> &dof_handler,
                                          const Vector<double> & solution) const
    {
        double point_value = 1e20;

        bool evaluation_point_found = false;
        for (const auto &cell : dof_handler.active_cell_iterators())
            if (!evaluation_point_found)
                for (const auto vertex : cell->vertex_indices())
                    if (cell->vertex(vertex).distance(evaluation_point) <
                        cell->diameter() * 1e-8)
                    {
                        point_value = solution(cell->vertex_dof_index(vertex, 0));

                        evaluation_point_found = true;
                        break;
                    }

        AssertThrow(evaluation_point_found,
                    ExcEvaluationPointNotFound(evaluation_point));

        std::cout << "   Point value=" << point_value << std::endl;
    }

    template <int dim>
    double
    PointValueEvaluation<dim>::operator()(const DoFHandler<dim> &dof_handler,
                                          const Vector<double> & solution)
    {
        double point_value = 1e20;

        bool evaluation_point_found = false;
        for (const auto &cell : dof_handler.active_cell_iterators())
            if (!evaluation_point_found)
                for (const auto vertex : cell->vertex_indices())
                    if (cell->vertex(vertex).distance(evaluation_point) <
                        cell->diameter() * 0.5)
                    {
                        point_value = solution(cell->vertex_dof_index(vertex, 0));

                        evaluation_point_found = true;
                        break;
                    }

        AssertThrow(evaluation_point_found,
                    ExcEvaluationPointNotFound(evaluation_point));

        return point_value;
    }

} // end namespace Evaluation


#endif //ION_PROPULSION_EVALUATION_H
