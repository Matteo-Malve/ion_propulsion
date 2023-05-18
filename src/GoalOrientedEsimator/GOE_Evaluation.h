#ifndef GETPOT_GOE_EVALUATION_H
#define GETPOT_GOE_EVALUATION_H

#include "GoalOrientedEstimator.h"

namespace GOE{
    using namespace dealii;
    namespace Evaluation{

        template <int dim>
        class EvaluationBase
        {
        public:
            virtual ~EvaluationBase() = default;

            void set_refinement_cycle(const unsigned int refinement_cycle);

            virtual void operator()(const DoFHandler<dim> &dof_handler,
                                    const Vector<double> & solution) const = 0;

        protected:
            unsigned int refinement_cycle;
        };



        template <int dim>
        void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step)
        {
            refinement_cycle = step;
        }


        template <int dim>
        class PointValueEvaluation : public EvaluationBase<dim>
        {
        public:
            PointValueEvaluation(const Point<dim> &evaluation_point);

            virtual void operator()(const DoFHandler<dim> &dof_handler,
                                    const Vector<double> & solution) const override;

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
        class PointXDerivativeEvaluation : public EvaluationBase<dim>
        {
        public:
            PointXDerivativeEvaluation(const Point<dim> &evaluation_point);

            virtual void operator()(const DoFHandler<dim> &dof_handler,
                                    const Vector<double> & solution) const;

            DeclException1(
                    ExcEvaluationPointNotFound,
                    Point<dim>,
            << "The evaluation point " << arg1
            << " was not found among the vertices of the present grid.");

        private:
            const Point<dim> evaluation_point;
        };


        template <int dim>
        PointXDerivativeEvaluation<dim>::PointXDerivativeEvaluation(
                const Point<dim> &evaluation_point)
                : evaluation_point(evaluation_point)
        {}


        template <int dim>
        void PointXDerivativeEvaluation<dim>::operator()(
                const DoFHandler<dim> &dof_handler,
                const Vector<double> & solution) const
        {
            double point_derivative = 0;

            QTrapezoid<dim>             vertex_quadrature;
            FEValues<dim>               fe_values(dof_handler.get_fe(),
                                                  vertex_quadrature,
                                                  update_gradients | update_quadrature_points);
            std::vector<Tensor<1, dim>> solution_gradients(vertex_quadrature.size());

            unsigned int evaluation_point_hits = 0;
            for (const auto &cell : dof_handler.active_cell_iterators())
                for (const auto vertex : cell->vertex_indices())
                    if (cell->vertex(vertex) == evaluation_point)
                    {
                        fe_values.reinit(cell);
                        fe_values.get_function_gradients(solution, solution_gradients);

                        unsigned int q_point = 0;
                        for (; q_point < solution_gradients.size(); ++q_point)
                            if (fe_values.quadrature_point(q_point) == evaluation_point)
                                break;

                        Assert(q_point < solution_gradients.size(), ExcInternalError());
                        point_derivative += solution_gradients[q_point][0];
                        ++evaluation_point_hits;

                        break;
                    }

            AssertThrow(evaluation_point_hits > 0,
                        ExcEvaluationPointNotFound(evaluation_point));

            point_derivative /= evaluation_point_hits;
            std::cout << "   Point x-derivative=" << point_derivative << std::endl;
        }




        template <int dim>
        class GridOutput : public EvaluationBase<dim>
        {
        public:
            GridOutput(const std::string &output_name_base);

            virtual void operator()(const DoFHandler<dim> &dof_handler,
                                    const Vector<double> & solution) const override;

        private:
            const std::string output_name_base;
        };


        template <int dim>
        GridOutput<dim>::GridOutput(const std::string &output_name_base)
                : output_name_base(output_name_base)
        {}


        template <int dim>
        void GridOutput<dim>::operator()(const DoFHandler<dim> &dof_handler,
                                         const Vector<double> & /*solution*/) const
        {
            std::ofstream out(output_name_base + "-" +
                              std::to_string(this->refinement_cycle) + ".svg");
            GridOut().write_svg(dof_handler.get_triangulation(), out);
        }



    } // end namespace GOE
} // end namespace Evaluation


#endif //GETPOT_GOE_EVALUATION_H