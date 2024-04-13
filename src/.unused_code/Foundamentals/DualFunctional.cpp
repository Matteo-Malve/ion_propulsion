#include "DualFunctional.h"

template <int dim>
void
EmitterFlux<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
                               Vector<double> &       rhs) const
{
    // Retrieve geometrical data from datafile
    double radius = redefined_3_datafile("wire_radius",250e-6);
    const Point<2> mesh_center(redefined_3_datafile("wire_center_x_coord",0.0),redefined_3_datafile("wire_center_y_coord",0.0));

    // Set up:
    rhs.reinit(dof_handler.n_dofs());
    // Quadrature
    QGauss<dim>          quadrature(dof_handler.get_fe().degree + 1);
    QGauss<dim-1>        face_quadrature(dof_handler.get_fe().degree + 1);
    // Finite elements
    FEValues<dim>        fe_values(dof_handler.get_fe(),
                                   quadrature,
                                   update_gradients | update_quadrature_points |
                                   update_JxW_values);
    FEFaceValues<dim>    fe_face_values(dof_handler.get_fe(),
                                        face_quadrature,
                                        update_gradients | update_quadrature_points |
                                        update_JxW_values | update_normal_vectors);

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
                if (mesh_center.distance(face_center) < 1.0000001 * radius) {
                    fe_face_values.reinit(cell, face);
                    cell->set_material_id(15);
                }
            }
    }

    // Now loop only on the previously marked cells
    for (const auto &cell : dof_handler.active_cell_iterators()) {
        fe_values.reinit(cell);
        if (cell->material_id() == 15) {
            cell_rhs = 0;
            fe_values.reinit(cell);
            // Retrieve the components of the normal vector to the boundary face
            Tensor<1, dim> n;
            for (const auto &face: cell->face_iterators())
                if (face->at_boundary()) {
                    const Point <dim> &face_center = face->center();
                    if (mesh_center.distance(face_center) < 1.0000001 * radius) {
                        fe_face_values.reinit(cell, face);
                        n = fe_face_values.normal_vector(0);
                    }
                }
            // Compute the flux for this cell
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
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                rhs(local_dof_indices[i]) += cell_rhs[i];
        }
    }
}

// #######################################
// Template initialization
// #######################################
template class EmitterFlux<2>;

