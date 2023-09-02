#ifndef GETPOT_DUALFUNCTIONAL_H
#define GETPOT_DUALFUNCTIONAL_H

#include "../includes&parameters_setup.h"

static GetPot redefined_3_datafile("../data_setup");

// The heart of this class. We inherit form a deal.ii class called Subscriptor
// It is a Base class containing tools for handling smart pointers
template <int dim>
class DualFunctionalBase : public Subscriptor
{
public:
    virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                              Vector<double> &       rhs) const = 0;
};


// Here we build the class to use our Goal Functional as rhs.
// Evidently, we could have skipped the double inheritance, but this way we open
// ourselves to the possibility to add new Goal Functionals in the future
template <int dim>
class EmitterFlux : public DualFunctionalBase<dim>
{
public:
    // Constructor
    EmitterFlux()=default;
    // The CORE
    virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                              Vector<double> &       rhs) const override;

};


template <int dim>
void
EmitterFlux<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
                               Vector<double> &       rhs) const
{
        double radius = redefined_3_datafile("Mesh/wire_radius",0.025);
        const Point<2> mesh_center(0, 0);
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
        const unsigned int dofs_per_face = dof_handler.get_fe().n_dofs_per_face();

        Vector<double>     face_rhs(dofs_per_face);
        std::vector<unsigned int> local_dof_face_indices(dofs_per_face);

        for (const auto &cell : dof_handler.active_cell_iterators()){
            fe_values.reinit(cell);
            for (const auto &face : cell->face_iterators())
                if (face->at_boundary()){
                    const Point<dim> & face_center = face->center();
                    if (mesh_center.distance(face_center)<1.00001*radius) {
                        fe_face_values.reinit(cell, face);
                        face_rhs = 0;

                        for (unsigned int q = 0; q < n_face_q_points; ++q) {
                            auto n = fe_face_values.normal_vector(q);
                            for (unsigned int i = 0; i < dofs_per_face; ++i) {
                                for (unsigned int k = 0; k < dim; ++k)
                                    face_rhs[i] += fe_face_values.shape_grad(i, q)[k] * (-n[k]);
                                face_rhs[i] *= fe_face_values.JxW(q);
                            }
                        }

                        face->get_dof_indices(local_dof_face_indices);
                        for (unsigned int i = 0; i < dofs_per_face; ++i)
                            rhs(local_dof_face_indices[i]) += face_rhs[i];
                    }
                }
        }

        //AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));

    }

#endif //GETPOT_DUALFUNCTIONAL_H
