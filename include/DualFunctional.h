#ifndef DUALFUNCTIONAL_H
#define DUALFUNCTIONAL_H

#include "includes.h"

namespace IonPropulsion{

	using namespace dealii;

  namespace DualFunctional
  {
    // @sect4{The DualFunctionalBase class}

    // First start with the base class for dual functionals. Since for linear
    // problems the characteristics of the dual problem play a role only in
    // the right hand side, we only need to provide for a function that
    // assembles the right hand side for a given discretization:
    template <int dim>
    class DualFunctionalBase : public Subscriptor
    {
    public:
      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
                                Vector<double> &       rhs) const = 0;
    };


    // @sect4{The dual functional PointValueEvaluation class}

    // As a first application, we consider the functional corresponding to the
    // evaluation of the solution's value at a given point which again we
    // assume to be a vertex. Apart from the constructor that takes and stores
    // the evaluation point, this class consists only of the function that
    // implements assembling the right hand side.
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


    // As for doing the main purpose of the class, assembling the right hand
    // side, let us first consider what is necessary: The right hand side of
    // the dual problem is a vector of values J(phi_i), where J is the error
    // functional, and phi_i is the i-th shape function. Here, J is the
    // evaluation at the point x0, i.e. J(phi_i)=phi_i(x0).
    //
    // Now, we have assumed that the evaluation point is a vertex. Thus, for
    // the usual finite elements we might be using in this program, we can
    // take for granted that at such a point exactly one shape function is
    // nonzero, and in particular has the value one. Thus, we set the right
    // hand side vector to all-zeros, then seek for the shape function
    // associated with that point and set the corresponding value of the right
    // hand side vector to one:
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


    // @sect4{The dual functional PointXDerivativeEvaluation class}

    // As second application, we again consider the evaluation of the
    // x-derivative of the solution at one point. Again, the declaration of
    // the class, and the implementation of its constructor is not too
    // interesting:
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


    // What is interesting is the implementation of this functional: here,
    // J(phi_i)=d/dx phi_i(x0).
    //
    // We could, as in the implementation of the respective evaluation object
    // take the average of the gradients of each shape function phi_i at this
    // evaluation point. However, we take a slightly different approach: we
    // simply take the average over all cells that surround this point. The
    // question which cells <code>surrounds</code> the evaluation point is
    // made dependent on the mesh width by including those cells for which the
    // distance of the cell's midpoint to the evaluation point is less than
    // the cell's diameter.
    //
    // Taking the average of the gradient over the area/volume of these cells
    // leads to a dual solution which is very close to the one which would
    // result from the point evaluation of the gradient. It is simple to
    // justify theoretically that this does not change the method
    // significantly.
    template <int dim>
    void PointXDerivativeEvaluation<dim>::assemble_rhs(
      const DoFHandler<dim> &dof_handler,
      Vector<double> &       rhs) const
    {
      // Again, first set all entries to zero:
      rhs.reinit(dof_handler.n_dofs());

      // Initialize a <code>FEValues</code> object with a quadrature formula,
      // have abbreviations for the number of quadrature points and shape
      // functions...
      QGauss<dim>        quadrature(dof_handler.get_fe().degree + 1);
      FEValues<dim>      fe_values(dof_handler.get_fe(),
                              quadrature,
                              update_gradients | update_quadrature_points |
                                update_JxW_values);
      const unsigned int n_q_points    = fe_values.n_quadrature_points;
      const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

      // ...and have two objects that are used to store the global indices of
      // the degrees of freedom on a cell, and the values of the gradients of
      // the shape functions at the quadrature points:
      Vector<double>            cell_rhs(dofs_per_cell);
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);

      // Finally have a variable in which we will sum up the area/volume of
      // the cells over which we integrate, by integrating the unit functions
      // on these cells:
      double total_volume = 0;

      // Then start the loop over all cells, and select those cells which are
      // close enough to the evaluation point:
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->center().distance(evaluation_point) <= cell->diameter())
          {
            // If we have found such a cell, then initialize the
            // <code>FEValues</code> object and integrate the x-component of
            // the gradient of each shape function, as well as the unit
            // function for the total area/volume.
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

      // After we have looped over all cells, check whether we have found any
      // at all, by making sure that their volume is non-zero. If not, then
      // the results will be botched, as the right hand side should then still
      // be zero, so throw an exception:
      AssertThrow(total_volume > 0,
                  ExcEvaluationPointNotFound(evaluation_point));

      // Finally, we have by now only integrated the gradients of the shape
      // functions, not taking their mean value. We fix this by dividing by
      // the measure of the volume over which we have integrated:
      rhs /= total_volume;
    }


  } // namespace DualFunctional

}
#endif