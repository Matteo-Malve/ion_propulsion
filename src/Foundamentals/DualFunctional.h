#ifndef GETPOT_DUALFUNCTIONAL_H
#define GETPOT_DUALFUNCTIONAL_H

#include "../includes&parameters_setup.h"

// ----------------------------------------------------------------------------


template <int dim>
class DualFunctionalBase : public Subscriptor
{
public:
    virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                              Vector<double> &       rhs) const = 0;
};

// ----------------------------------------------------------------------------

template <int dim>
class EmitterFlux : public DualFunctionalBase<dim>
{
public:
    EmitterFlux()=default;

    virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                              Vector<double> &       rhs) const override;
    /*void
    assemble_rhs_with_dealii_n(const DoFHandler<dim> &dof_handler,
                               Vector<double> &       rhs) const;*/

};

// ----------------------------------------------------------------------------

/*
template <int dim>
void
EmitterFlux<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
                               Vector<double> &       rhs) const
{
    double wire_radius = 0.025;     // Terribile, da generalizzare
    rhs.reinit(dof_handler.n_dofs());
    QGauss<dim>          quadrature(dof_handler.get_fe().degree + 1);
    QGauss<dim-1>        face_quadrature(dof_handler.get_fe().degree + 1);
    
    FEValues<dim>        fe_values(dof_handler.get_fe(),
                                 quadrature,
                                 update_gradients | update_quadrature_points |
                                 update_JxW_values);  

    FEFaceValues<dim>    fe_face_values(dof_handler.get_fe(),
                                        face_quadrature,
                                        update_gradients | update_quadrature_points |
                                        update_JxW_values | update_normal_vectors);

    const unsigned int n_face_q_points = face_quadrature.size();
    const unsigned int n_q_points      = quadrature.size();
    const unsigned int dofs_per_face = dof_handler.get_fe().n_dofs_per_face();
    const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();

    Vector<double>     cell_rhs_face(dofs_per_face);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<unsigned int> local_dof_face_indices(dofs_per_face);
    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
       for (const auto &face : cell->face_iterators())
          if (face->at_boundary()){
            const Point<dim> c = face->center();
            if ((c[1] < 0.2) && (c[0] <= wire_radius) && (c[0] >= -wire_radius))

                cell_rhs = 0;
            for (unsigned int q = 0; q < n_face_q_points; ++q) {
                for (const auto i : face->vertex_indices()) {
                    auto v = face->vertex(i);
                    double theta = std::atan2(v[1]/v[0]);
                    std::vector<double> n = {cos(theta), sin(theta)};
                    for (unsigned int k = 0; k < dim; ++k)
                        cell_rhs(i) += fe_face_values.shape_grad(i, q)[k] * n[k];
                    cell_rhs(i) = cell_rhs(i) * fe_face_values.JxW(q);
                }
            }
            face->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_face; ++i)
                rhs(local_dof_indices[i]) += cell_rhs(i);
        }

    //AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));

}

 */

template <int dim>
void
EmitterFlux<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
                               Vector<double> &       rhs) const
{
    double wire_radius = 0.025;     // Terribile, da generalizzare
    rhs.reinit(dof_handler.n_dofs());
    QGauss<dim>          quadrature(dof_handler.get_fe().degree + 1);
    QGauss<dim-1>        face_quadrature(dof_handler.get_fe().degree + 1);
    
    FEValues<dim>        fe_values(dof_handler.get_fe(),
                                 quadrature,
                                 update_gradients | update_quadrature_points |
                                 update_JxW_values);  

    FEFaceValues<dim>    fe_face_values(dof_handler.get_fe(),
                                        face_quadrature,
                                        update_gradients | update_quadrature_points |
                                        update_JxW_values | update_normal_vectors);

    const unsigned int n_face_q_points = face_quadrature.size();
    // const unsigned int n_q_points      = quadrature.size();
    const unsigned int dofs_per_face = dof_handler.get_fe().n_dofs_per_face();
    //const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();

    //Vector<double>     cell_rhs_face(dofs_per_face);
    Vector<double>     cell_rhs(dofs_per_face);
    std::vector<unsigned int> local_dof_face_indices(dofs_per_face);
    //std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        for (const auto &face : cell->face_iterators())
          if (face->at_boundary()){
            const Point<dim> c = face->center();
            if ((c[1] < 0.2) && (c[0] <= wire_radius) && (c[0] >= -wire_radius)){
                fe_face_values.reinit(cell, face);
                cell_rhs = 0;   // ?    Vettore = 0 ?
                                // Perchè mettiamo a zero in questo rettangolo?
            for (unsigned int q = 0; q < n_face_q_points; ++q) {
                auto n = fe_face_values.normal_vector(q);
                for (unsigned int i = 0; i < dofs_per_face; ++i) {
                    for (unsigned int k = 0; k < dim; ++k)
                        cell_rhs[i] += fe_face_values.shape_grad(i, q)[k] * n[k];
                    cell_rhs[i] *= fe_face_values.JxW(q);  // ??    access operator [], perchè () ?
                }
            }
          
            face->get_dof_indices(local_dof_face_indices);
            for (unsigned int i = 0; i < dofs_per_face; ++i)
                rhs(local_dof_face_indices[i]) += cell_rhs[i];   // ??  idem sopra
        }
    }
    }

    //AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));

}   

#endif //GETPOT_DUALFUNCTIONAL_H
