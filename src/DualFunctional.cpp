#include "DualFunctional.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace DualFunctional{
    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------
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
      // So, first set everything to zeros...
      rhs.reinit(dof_handler.n_dofs());

      // ...then loop over cells and find the evaluation point among the
      // vertices (or very close to a vertex, which may happen due to floating
      // point round-off):
      for (const auto &cell : dof_handler.active_cell_iterators())
        for (const auto vertex : cell->vertex_indices())
          if (cell->vertex(vertex).distance(evaluation_point) <
              cell->diameter() * 1e-8)
          {
            // Ok, found, so set corresponding entry, and leave function
            // since we are finished:
            rhs(cell->vertex_dof_index(vertex, 0)) = 1;
            return;
          }

      // Finally, a sanity check: if we somehow got here, then we must have
      // missed the evaluation point, so raise an exception unconditionally:
      AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));
    }

    // ------------------------------------------------------
    // PointXDerivativeEvaluation
    // ------------------------------------------------------
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

      // Then start the loop over all cells, and select those cells which are
      // close enough to the evaluation point:
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->center().distance(evaluation_point) <= cell->diameter())
          {

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

            // If we have the local contributions, distribute them to the
            // global vector:
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              rhs(local_dof_indices[i]) += cell_rhs(i);
          }
      AssertThrow(total_volume > 0, ExcEvaluationPointNotFound(evaluation_point));
      rhs /= total_volume;
    }

    // Template instantiation
    template class DualFunctionalBase<2>;
    template class PointValueEvaluation<2>;
    template class PointXDerivativeEvaluation<2>;

  } // namespace DualFunctional
} // namespace IonPropulsion