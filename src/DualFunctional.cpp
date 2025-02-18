#include "DualFunctional.h"
#include <Refinement.h>

namespace IonPropulsion{
  using namespace dealii;
  namespace DualFunctional{
    // ------------------------------------------------------
    // DualFunctionalBase
    // ------------------------------------------------------
    template <int dim>
    DualFunctionalBase<dim>::DualFunctionalBase(const unsigned mapping_degree):
      mapping(mapping_degree)
    {}

    // ------------------------------------------------------
    // PointValueEvaluation
    // ------------------------------------------------------
    template <int dim>
      PointValueEvaluation<dim>::PointValueEvaluation(
        const unsigned mapping_degree,
        const Point<dim> &evaluation_point)
        : DualFunctionalBase<dim>(mapping_degree), evaluation_point(evaluation_point)
    {}

    template <int dim>
    void
    PointValueEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &,
                                            PETScWrappers::MPI::Vector &
                                            , AffineConstraints<double> &) const
    { //TODO
      /*
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
          if (distance < cell->diameter() * 1.e-8)
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
      }*/
    }

    // ------------------------------------------------------
    // StandardFluxEvaluation
    // ------------------------------------------------------
    template <int dim>
      StandardFluxEvaluation<dim>::StandardFluxEvaluation(
        const unsigned mapping_degree,
        const std::set<unsigned int> &boundary_ids)
        : DualFunctionalBase<dim>(mapping_degree), boundary_ids(boundary_ids)
    {}

    template <int dim>
    void
    StandardFluxEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
                                            PETScWrappers::MPI::Vector &       rhs
                                            , AffineConstraints<double> & hanging_node_constraints) const
    {

      // Set up:
      //rhs.reinit(dof_handler.n_dofs());  // Done already by the calling function in solver

      auto & fe_face = dof_handler.get_fe();

      // Quadrature
      const QGauss<dim-1> face_quadrature(fe_face.degree + 1);

      // Finite elements
      FEFaceValues<dim> fe_face_values(this->mapping, fe_face,
                                       face_quadrature,
                                       update_gradients |
                                       update_normal_vectors |
                                       update_JxW_values);

      const unsigned int dofs_per_cell = fe_face.n_dofs_per_cell();
      const unsigned int n_face_q_points = face_quadrature.size();
      Vector<double> cell_rhs(dofs_per_cell);
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);

      // Loop over all cells and faces
      for (const auto &cell : dof_handler.active_cell_iterators()) {
        if (cell->is_locally_owned()) {
          cell_rhs = 0;

          for (const auto &face : cell->face_iterators()) {
            // Process only boundary faces with the specified boundary_id
            if (face->at_boundary() && boundary_ids.count(face->boundary_id()) > 0) {
              fe_face_values.reinit(cell, face);

              // Compute flux for this face
              for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) {
                const Tensor<1, dim> &n = fe_face_values.normal_vector(q_point);

                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                  // Sum the flux contribution from each quadrature point on the boundary
                  cell_rhs[i] += (fe_face_values.shape_grad(i, q_point) * (-n)) * fe_face_values.JxW(q_point);
                }
              }
            }
          }

          // Local to global: sum the contribution of this cell to the rhs
          cell->get_dof_indices(local_dof_indices);
          hanging_node_constraints.distribute_local_to_global(
            cell_rhs,
            local_dof_indices,
            rhs
            );
        }
      }
    }





    // Template instantiation
    template class DualFunctionalBase<2>;
    template class PointValueEvaluation<2>;
    template class StandardFluxEvaluation<2>;


  } // namespace DualFunctional
} // namespace IonPropulsion