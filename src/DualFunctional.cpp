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

      // Keep track of nearest vertex
      double min_distance = std::numeric_limits<double>::max();
      typename DoFHandler<dim>::active_cell_iterator nearest_cell;
      unsigned int nearest_vertex = 0;

      // ...then loop over cells and find the evaluation point among the
      // vertices (or very close to a vertex, which may happen due to floating
      // point round-off):
      for (const auto &cell : dof_handler.active_cell_iterators())
        for (const auto vertex : cell->vertex_indices()) {
          double distance = cell->vertex(vertex).distance(evaluation_point);
          if (distance < cell->diameter() * 1.)
          {
            // Ok, found, so set corresponding entry, and leave function
            // since we are finished:
            rhs(cell->vertex_dof_index(vertex, 0)) = 1;
            return;
          }
          // Keep track of the nearest vertex if the exact point isn't found
          if (distance < min_distance) {
            min_distance = distance;
            nearest_cell = cell;
            nearest_vertex = vertex;
          }
        }
      std::cout<<"      DualRhs: Point not found" << std::endl;
      // If no vertex within tolerance was found, set rhs at the nearest vertex
      if (min_distance < std::numeric_limits<double>::max()) {
        rhs(nearest_cell->vertex_dof_index(nearest_vertex, 0)) = 1;
      } else {
        AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));
      }
    }

    // ------------------------------------------------------
    // FluxEvaluation
    // ------------------------------------------------------
    template <int dim>
      FluxEvaluation<dim>::FluxEvaluation(
        const unsigned int boundary_id)
        : boundary_id(boundary_id)
    {}

    template <int dim>
    void
    FluxEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
                                            Vector<double> &       rhs) const
    {
      // Set up:
      rhs.reinit(dof_handler.n_dofs());

      auto & fe_face = dof_handler.get_fe();

      // Quadrature
      const QGauss<dim-1> face_quadrature(fe_face.degree + 1);

      // Finite elements
      FEFaceValues<dim> fe_face_values(fe_face,
                                       face_quadrature,
                                       update_gradients |
                                       update_normal_vectors |
                                       update_JxW_values);

      const unsigned int dofs_per_cell = fe_face.n_dofs_per_cell();
      const unsigned int n_face_q_points = face_quadrature.size();
      const unsigned int n_facet_dofs = fe_face.n_dofs_per_cell();
      Vector<double> cell_rhs(dofs_per_cell);
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);

      // Loop over all cells and faces
      for (const auto &cell : dof_handler.active_cell_iterators()) {
        cell_rhs = 0;

        for (const auto &face : cell->face_iterators()) {
          // Process only boundary faces with the specified boundary_id
          if (face->at_boundary() && face->boundary_id() == 1) {
            fe_face_values.reinit(cell, face);

            // Compute flux for this face
            for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) {
              const Tensor<1, dim> &n = fe_face_values.normal_vector(q_point);

              for (unsigned int i = 0; i < n_facet_dofs; ++i) {
                // Sum the flux contribution from each quadrature point on the boundary
                cell_rhs[i] += (fe_face_values.shape_grad(i, q_point) * (-n)) * fe_face_values.JxW(q_point);
              }
            }
          }
        }

        // Local to global: sum the contribution of this cell to the rhs
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < n_facet_dofs; ++i)
          rhs(local_dof_indices[i]) += cell_rhs[i];
      }
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
    template class FluxEvaluation<2>;
    template class PointXDerivativeEvaluation<2>;

  } // namespace DualFunctional
} // namespace IonPropulsion