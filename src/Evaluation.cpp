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

		// ------------------------------------------------------
		// PointValueEvaluation
		// ------------------------------------------------------

		template <int dim>
		PointValueEvaluation<dim>::PointValueEvaluation(
			const Point<dim> &evaluation_point)
			: evaluation_point(evaluation_point)
		{}

		template <int dim>
		std::pair<std::string, double>
		PointValueEvaluation<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const Vector<double> & solution,
			const Triangulation<dim> &       triangulation) const
		{
			(void)triangulation;
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

			//AssertThrow(evaluation_point_found, ExcEvaluationPointNotFound(evaluation_point));
			if (!evaluation_point_found) {
				cout<<"        Evaluation: Vertex not found in the mesh."<<std::endl;
				point_value = VectorTools::point_value(dof_handler, solution, evaluation_point);
			}

			std::cout << std::scientific << std::setprecision(12)
								<< "   Point value = " << point_value << std::endl;
			std::cout << std::scientific << std::setprecision(12)
								<< "   Exact value = " << EXACT_POINT_VALUE << std::endl;

			// Update table with exact error
			double exact_error = std::fabs(point_value-EXACT_POINT_VALUE);
			return std::make_pair("ex POINT err",exact_error);
		}

		// ------------------------------------------------------
		// FluxEvaluation
		// ------------------------------------------------------

		template <int dim>
		FluxEvaluation<dim>::FluxEvaluation()
		{}

		template <int dim>
		std::pair<std::string, double>
		FluxEvaluation<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const Vector<double>  &solution,
			const Triangulation<dim> &       triangulation) const
		{
			(void)triangulation;
			double flux = 0;
			const QGauss<dim-1> face_quadrature(7);			// dof_handler.get_fe().degree + 1
			FEFaceValues<dim> fe_face_values(dof_handler.get_fe(),
																			 face_quadrature,
																			 update_gradients | update_normal_vectors | update_JxW_values);
			const unsigned int n_face_q_points = face_quadrature.size();
			std::vector<Tensor<1, dim>> solution_gradients(n_face_q_points);

			for (const auto &cell : dof_handler.active_cell_iterators())
				for (const auto &face : cell->face_iterators())
					if (face->at_boundary() && face->boundary_id() == 1) {
						fe_face_values.reinit(cell, face);
						fe_face_values.get_function_gradients(solution, solution_gradients);
						for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) {
							const Tensor<1, dim> &n = fe_face_values.normal_vector(q_point);
							flux += ((-solution_gradients[q_point]) * (-n)) * fe_face_values.JxW(q_point);
						}
					}

			std::cout << std::scientific << std::setprecision(12)
								<< "   Comp. Flux = " << flux << std::endl;
			std::cout << std::scientific << std::setprecision(12)
								<< "   Exact Flux = " << EXACT_FLUX << std::endl;

			// Update table with exact error
			double exact_error = std::fabs(flux-EXACT_FLUX);
			return std::make_pair("std FLUX err",exact_error);
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
		GridOutput<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const Vector<double> & /*solution*/,
			const Triangulation<dim> &       triangulation) const
		{
			(void)triangulation;
			std::ofstream out(output_name_base + "-" +
												std::to_string(this->refinement_cycle) + ".svg");
			GridOut().write_svg(dof_handler.get_triangulation(), out);
			return std::make_pair("null",-1.0e-20);
		}

		// ------------------------------------------------------
		// L2_error_estimate
		// ------------------------------------------------------

		template <int dim>
		L2_error_estimate<dim>::L2_error_estimate(const Function<dim> & analytical_solution)
			: analytical_solution(&analytical_solution)
		{}


		template <int dim>
		std::pair<std::string, double>
		L2_error_estimate<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const Vector<double> & solution,
			const Triangulation<dim> &       triangulation) const
		{
			// Set up an error vector
			Vector<double> difference_per_cell(triangulation.n_active_cells());

			// Compute the difference between the finite element solution and the exact solution
			const QGauss<dim> quadrature_formula(2*dof_handler.get_fe().degree + 1);
			VectorTools::integrate_difference(
					dof_handler,
					solution,
					*analytical_solution,
					difference_per_cell,
					quadrature_formula,
					VectorTools::L2_norm);

			const double L2_error = VectorTools::compute_global_error(triangulation,
																																difference_per_cell,
																																VectorTools::L2_norm);
			return std::make_pair("L2",L2_error);
		}

		// ------------------------------------------------------
		// H1_error_estimate
		// ------------------------------------------------------

		template <int dim>
		H1_error_estimate<dim>::H1_error_estimate(const Function<dim> & analytical_solution)
			: analytical_solution(&analytical_solution)
		{}


		template <int dim>
		std::pair<std::string, double>
		H1_error_estimate<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const Vector<double> & solution,
			const Triangulation<dim> &       triangulation) const
		{
			// Set up an error vector
			Vector<double> difference_per_cell(triangulation.n_active_cells());

			// Compute the difference between the finite element solution and the exact solution
			const QGauss<dim> quadrature_formula(2*dof_handler.get_fe().degree + 1);
			VectorTools::integrate_difference(
					dof_handler,
					solution,
					*analytical_solution,
					difference_per_cell,
					quadrature_formula,
					VectorTools::H1_norm);

			const double H1_error = VectorTools::compute_global_error(triangulation,
																																difference_per_cell,
																																VectorTools::H1_norm);
			return std::make_pair("H1",H1_error);
		}


		// Template instantiation
		template class EvaluationBase<2>;
		template class PointValueEvaluation<2>;
		template class FluxEvaluation<2>;
		template class L2_error_estimate<2>;
		template class H1_error_estimate<2>;
		template class GridOutput<2>;


	} // namespace Evaluation
} // namespace IonPropulsion