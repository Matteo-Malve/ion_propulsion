#include "Evaluation.h"

#include <Constants.h>

namespace IonPropulsion{
	using namespace dealii;
	namespace Evaluation{
		// ------------------------------------------------------
		// EvaluationBase
		// ------------------------------------------------------
		template <int dim>
		void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step)
		{
			refinement_cycle = step;
		}

		template <int dim>
		PointValueEvaluation<dim>::PointValueEvaluation(
			const Point<dim> &evaluation_point)
			: evaluation_point(evaluation_point)
		{}

		// ------------------------------------------------------
		// PointValueEvaluation
		// ------------------------------------------------------

		template <int dim>
		std::pair<std::string, double>
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

			std::cout << std::scientific << std::setprecision(12)
								<< "   Point value=" << point_value << std::endl;

			// Update table with exact error
			double exact_error = std::fabs(point_value-EXACT_POINT_VALUE);
			return std::make_pair("ex POINT err",exact_error);
		}

		// ------------------------------------------------------
		// PointXDerivativeEvaluation
		// ------------------------------------------------------

		template <int dim>
    PointXDerivativeEvaluation<dim>::PointXDerivativeEvaluation(
      const Point<dim> &evaluation_point)
      : evaluation_point(evaluation_point)
    {}


    // The more interesting things happen inside the function doing the actual
    // evaluation:
    template <int dim>
		std::pair<std::string, double>
		PointXDerivativeEvaluation<dim>::operator()(
      const DoFHandler<dim> &dof_handler,
      const Vector<double> & solution) const
    {
      // This time initialize the return value with something useful, since we
      // will have to add up a number of contributions and take the mean value
      // afterwards...
      double point_derivative = 0;

      // ...then have some objects of which the meaning will become clear
      // below...
      QTrapezoid<dim>             vertex_quadrature;
      FEValues<dim>               fe_values(dof_handler.get_fe(),
                              vertex_quadrature,
                              update_gradients | update_quadrature_points);
      std::vector<Tensor<1, dim>> solution_gradients(vertex_quadrature.size());

      // ...and next loop over all cells and their vertices, and count how
      // often the vertex has been found:
      unsigned int evaluation_point_hits = 0;
      for (const auto &cell : dof_handler.active_cell_iterators())
        for (const auto vertex : cell->vertex_indices())
          if (cell->vertex(vertex) == evaluation_point)
            {
              // Things are now no more as simple, since we can't get the
              // gradient of the finite element field as before, where we
              // simply had to pick one degree of freedom at a vertex.
              //
              // Rather, we have to evaluate the finite element field on this
              // cell, and at a certain point. As you know, evaluating finite
              // element fields at certain points is done through the
              // <code>FEValues</code> class, so we use that. The question is:
              // the <code>FEValues</code> object needs to be a given a
              // quadrature formula and can then compute the values of finite
              // element quantities at the quadrature points. Here, we don't
              // want to do quadrature, we simply want to specify some points!
              //
              // Nevertheless, the same way is chosen: use a special
              // quadrature rule with points at the vertices, since these are
              // what we are interested in. The appropriate rule is the
              // trapezoidal rule, so that is the reason why we used that one
              // above.
              //
              // Thus: initialize the <code>FEValues</code> object on this
              // cell,
              fe_values.reinit(cell);
              // and extract the gradients of the solution vector at the
              // vertices:
              fe_values.get_function_gradients(solution, solution_gradients);

              // Now we have the gradients at all vertices, so pick out that
              // one which belongs to the evaluation point (note that the
              // order of vertices is not necessarily the same as that of the
              // quadrature points):
              unsigned int q_point = 0;
              for (; q_point < solution_gradients.size(); ++q_point)
                if (fe_values.quadrature_point(q_point) == evaluation_point)
                  break;

              // Check that the evaluation point was indeed found,
              Assert(q_point < solution_gradients.size(), ExcInternalError());
              // and if so take the x-derivative of the gradient there as the
              // value which we are interested in, and increase the counter
              // indicating how often we have added to that variable:
              point_derivative += solution_gradients[q_point][0];
              ++evaluation_point_hits;

              // Finally break out of the innermost loop iterating over the
              // vertices of the present cell, since if we have found the
              // evaluation point at one vertex it cannot be at a following
              // vertex as well:
              break;
            }

      // Now we have looped over all cells and vertices, so check whether the
      // point was found:
      AssertThrow(evaluation_point_hits > 0,
                  ExcEvaluationPointNotFound(evaluation_point));

      // We have simply summed up the contributions of all adjacent cells, so
      // we still have to compute the mean value. Once this is done, report
      // the status:
      point_derivative /= evaluation_point_hits;
      std::cout << "   Point x-derivative=" << point_derivative << std::endl;
			return std::make_pair("null",-1.0e-20);
    }

		// ------------------------------------------------------
		// GridOutput
		// ------------------------------------------------------

		template <int dim>
		GridOutput<dim>::GridOutput(const std::string &output_name_base)
			: output_name_base(output_name_base)
		{}


		template <int dim>
		std::pair<std::string, double>
		GridOutput<dim>::operator()(const DoFHandler<dim> &dof_handler,
																		 const Vector<double> & /*solution*/) const
		{
			std::ofstream out(output_name_base + "-" +
												std::to_string(this->refinement_cycle) + ".svg");
			GridOut().write_svg(dof_handler.get_triangulation(), out);
			return std::make_pair("null",-1.0e-20);
		}


		// Template instantiation
		template class GridOutput<2>;
		template class PointXDerivativeEvaluation<2>;
		template class PointValueEvaluation<2>;
		template class EvaluationBase<2>;

	} // namespace Evaluation
} // namespace IonPropulsion