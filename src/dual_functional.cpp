#include "dual_functional.h"


// ------------------------------------------------------------      
// POINT evaluation
// ------------------------------------------------------------      

template <int dim>
PointValueEvaluation<dim>::PointValueEvaluation(const Point<dim> &evaluation_point): evaluation_point(evaluation_point)
{}

template <int dim>
void PointValueEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,Vector<double> &rhs) const{
  rhs.reinit(dof_handler.n_dofs());

  // Keep track of nearest vertex
  double min_distance = std::numeric_limits<double>::max();
  typename DoFHandler<dim>::active_cell_iterator nearest_cell;
  unsigned int nearest_vertex = 0;

  for(const auto &cell : dof_handler.active_cell_iterators())
    for(const auto vertex : cell->vertex_indices()){
      double distance = cell->vertex(vertex).distance(evaluation_point);

      // Check if the vertex is close enough to the evaluation point
      if (distance < cell->diameter() * 1e-8) {
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
    
  // If no vertex within tolerance was found, set rhs at the nearest vertex
  if (min_distance < std::numeric_limits<double>::max()) {
    rhs(nearest_cell->vertex_dof_index(nearest_vertex, 0)) = 1;
  } else {
    AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));
  }
}

// ------------------------------------------------------------      
// POINT dY
// ------------------------------------------------------------      

template <int dim>
PointYDerivativeEvaluation<dim>::PointYDerivativeEvaluation(
  const Point<dim> &evaluation_point)
  : evaluation_point(evaluation_point)
{}

template <int dim>
void PointYDerivativeEvaluation<dim>::assemble_rhs(
  const DoFHandler<dim> &dof_handler,
  Vector<double>        &rhs) const
{
  rhs.reinit(dof_handler.n_dofs());
  const QGauss<dim>  quadrature(dof_handler.get_fe().degree + 1);
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
      {
        fe_values.reinit(cell);
        cell_rhs = 0;
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_rhs(i) +=
                fe_values.shape_grad(i, q)[1] // (d/dy phi_i(x_q))
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

// ------------------------------------------------------------      
// AREA evaluation
// ------------------------------------------------------------      


template <int dim>
AreaEvaluation<dim>::AreaEvaluation(
  const Point<dim> &center_point,
  const double radius)
  : center_point(center_point), radius(radius)
{}

template <int dim>
void AreaEvaluation<dim>::assemble_rhs(
  const DoFHandler<dim> &dof_handler,
  Vector<double>        &rhs) const
{
  rhs.reinit(dof_handler.n_dofs());
  const QGauss<dim>  quadrature(dof_handler.get_fe().degree + 1);
  FEValues<dim>      fe_values(dof_handler.get_fe(),
                          quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int n_q_points    = fe_values.n_quadrature_points;
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  Vector<double>            cell_rhs(dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  double total_volume = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->center().distance(center_point) <= radius)
      {
        fe_values.reinit(cell);
        cell_rhs = 0;
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_rhs(i) +=
                fe_values.shape_value(i, q) // phi_i(x_q)
                * fe_values.JxW(q);           // * dx
            total_volume += fe_values.JxW(q);
          }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          rhs(local_dof_indices[i]) += cell_rhs(i);
      }
  AssertThrow(total_volume > 0,
              ExcEvaluationPointNotFound(center_point));
  rhs /= total_volume;
}

// ------------------------------------------------------------      
// BOUNDARY FLUX
// ------------------------------------------------------------      

template <int dim>
BoundaryFluxEvaluation<dim>::BoundaryFluxEvaluation(const unsigned int boundary_id)
  : boundary_id(boundary_id)
{}

template <int dim>
void BoundaryFluxEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler, Vector<double> &rhs) const {
  // Set up:
  rhs.reinit(dof_handler.n_dofs());
  // Quadrature
  const QGauss<dim>          quadrature(dof_handler.get_fe().degree + 1);
  const QGauss<dim-1>        face_quadrature(dof_handler.get_fe().degree + 1);  
  // Finite elements
  FEValues<dim>       fe_values(dof_handler.get_fe(),
                                quadrature,
                                update_gradients | update_quadrature_points | update_JxW_values);
  FEFaceValues<dim>   fe_face_values(dof_handler.get_fe(), 
                                     face_quadrature,
                                     update_gradients | update_normal_vectors | update_JxW_values);
  /*
  const unsigned int dofs_per_cell = fe_face_values.get_fe().n_dofs_per_cell();
  const unsigned int n_face_q_points = face_quadrature.size();
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  // Loop over all cells
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    cell_rhs = 0;

    // Loop over all faces of the cell
    for (const auto &face : cell->face_iterators()) {
      if (face->at_boundary() && face->boundary_id() == boundary_id) {
        fe_face_values.reinit(cell, face);

        // Compute flux for this cell
        for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) {
          for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            cell_rhs(i) += fe_face_values.shape_grad(i, q_point) *      // grad_phi_i
                           (-fe_face_values.normal_vector(q_point)) *   // - normal_vector (inwards)
                           fe_face_values.JxW(q_point);                 // d(gamma)
          }
        }
      }
    }

    // Distribute local to global
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }*/
  
  const unsigned int n_q_points = quadrature.size();
  const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // Loop over all cells and mark those belonging to the first foil around the emitter
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    fe_values.reinit(cell);
    for (const auto &face: cell->face_iterators())
      if (face->at_boundary()) {
        const Point <dim> &face_center = face->center();
        if (face_center[1] < 1.0000001 * l   &&    std::abs(face_center[0]) < 1.0000001 * l) {
          fe_face_values.reinit(cell, face);
          cell->set_material_id(15);
        }
      }
  }

  // Now loop only on the previously marked cells
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    bool flag = false;
    // Check if cell is of interest
    for (const auto &face: cell->face_iterators())
      if (face->at_boundary()    &&    face->center()[1] < 1.0000001 * l   &&    std::abs(face->center()[0]) < 1.0000001 * l)
        flag = true;
    // Cell passed the test
    if (flag) {
      cell_rhs = 0;
      fe_values.reinit(cell);
      // [FEFaceValues] Retrieve the components of the normal vector to the boundary face
      Tensor<1, dim> n;
      for (const auto &face: cell->face_iterators())
        if (face->at_boundary())
          if (face->at_boundary()    &&    face->center()[1] < 1.0000001 * l   &&    std::abs(face->center()[0]) < 1.0000001 * l) {
              fe_face_values.reinit(cell, face);
              n = fe_face_values.normal_vector(0);
            }    
      // [FEValues] Compute the flux for this cell
      for (unsigned int q = 0; q < n_q_points; ++q) {
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
          for (unsigned int k = 0; k < dim; ++k) {
              cell_rhs[i] += fe_values.shape_grad(i, q)[k] * (-n[k]);
          }
          cell_rhs[i] *= fe_values.JxW(q);
        }
      }
      // Local to Global: sum the contribution of this cell to the rhs
      cell->get_dof_indices(local_dof_indices);
      rhs.add(local_dof_indices,cell_rhs);
    }
  }
}


template <int dim>
FaceBoundaryFluxEvaluation<dim>::FaceBoundaryFluxEvaluation(const unsigned int boundary_id)
  : boundary_id(boundary_id)
{}

template <int dim>
void FaceBoundaryFluxEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler, Vector<double> &rhs) const {
  // Set up:
  rhs.reinit(dof_handler.n_dofs());

  // Quadrature
  const QGauss<dim-1> face_quadrature(dof_handler.get_fe().degree + 1);  

  // Finite elements
  FEFaceValues<dim> fe_face_values(dof_handler.get_fe(), 
                                   face_quadrature,
                                   update_gradients | update_normal_vectors | update_JxW_values);

  const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();
  const unsigned int n_face_q_points = face_quadrature.size();
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
          
          for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            // Sum the flux contribution from each quadrature point on the boundary
            cell_rhs[i] += (fe_face_values.shape_grad(i, q_point) * (-n)) * fe_face_values.JxW(q_point);
          }
        }
      }
    }

    // Local to global: sum the contribution of this cell to the rhs
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      rhs(local_dof_indices[i]) += cell_rhs[i];
  }
}


// #######################################
// Template initialization
// #######################################
template class PointValueEvaluation<2>;
template class PointYDerivativeEvaluation<2>;
template class AreaEvaluation<2>;
template class BoundaryFluxEvaluation<2>;
template class FaceBoundaryFluxEvaluation<2>;

