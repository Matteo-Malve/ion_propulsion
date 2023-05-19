
#ifndef GETPOT_GOE_DUALFUNCTIONAL_H
#define GETPOT_GOE_DUALFUNCTIONAL_H

#include "GOE_Data.h"

namespace GOE{
    using namespace dealii;
    namespace DualFunctional
    {

        template <int dim>
        class DualFunctionalBase : public Subscriptor
        {
        public:
            virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                      Vector<double> &       rhs) const = 0;
        };



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


        template <int dim>
        PointValueEvaluation<dim>::PointValueEvaluation(
                const Point<dim> &evaluation_point)
                : evaluation_point(evaluation_point)
        {}


        template <int dim>
        void
        PointValueEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
                                                Vector<double> &       rhs) const
        {
            rhs.reinit(dof_handler.n_dofs());

            for (const auto &cell : dof_handler.active_cell_iterators())
                for (const auto vertex : cell->vertex_indices())
                    if (cell->vertex(vertex).distance(evaluation_point) <
                        cell->diameter() * 1e-8)
                    {
                        rhs(cell->vertex_dof_index(vertex, 0)) = 1;
                        return;
                    }

            AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));
        }



        template <int dim>
        class PointXDerivativeEvaluation : public DualFunctionalBase<dim>
        {
        public:
            PointXDerivativeEvaluation(const Point<dim> &evaluation_point);

            virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                      Vector<double> &       rhs) const;

            DeclException1(
                    ExcEvaluationPointNotFound,
                    Point<dim>,
            << "The evaluation point " << arg1
            << " was not found among the vertices of the present grid.");

        protected:
            const Point<dim> evaluation_point;
        };


        template <int dim>
        PointXDerivativeEvaluation<dim>::PointXDerivativeEvaluation(
                const Point<dim> &evaluation_point)
                : evaluation_point(evaluation_point)
        {}


        template <int dim>
        void PointXDerivativeEvaluation<dim>::assemble_rhs(
                const DoFHandler<dim> &dof_handler,
                Vector<double> &       rhs) const
        {
            rhs.reinit(dof_handler.n_dofs());

            QGauss<dim>        quadrature(dof_handler.get_fe().degree + 1);
            FEValues<dim>      fe_values(dof_handler.get_fe(),
                                         quadrature,
                                         update_gradients | update_quadrature_points |
                                         update_JxW_values);
            const unsigned int n_q_points    = fe_values.n_quadrature_points;
            const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

            Vector<double>            cell_rhs(dofs_per_cell);
            std::vector<unsigned int> local_dof_indices(dofs_per_cell);

            double total_volume = 0;

            for (const auto &cell : dof_handler.active_cell_iterators())
                if (cell->center().distance(evaluation_point) <= cell->diameter())
                {fcdxc
                    fe_values.reinit(cell);
                    cell_rhs = 0;

                    for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                            cell_rhs(i) +=
                                    fe_values.shape_grad(i, q)[0] // (d/dx phi_i(x_q))
                                    * fe_values.JxW(q);           // * dx
                        total_volume += fe_values.JxW(q);
                    }

                    cell->get_dof_indices(local_dof_indices);
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        rhs(local_dof_indices[i]) += cell_rhs(i);
                }

            AssertThrow(total_volume > 0,
                        ExcEvaluationPointNotFound(evaluation_point));

            rhs /= total_volume;
        }


    } // namespace DualFunctional
}

#endif //GETPOT_GOE_DUALFUNCTIONAL_H
