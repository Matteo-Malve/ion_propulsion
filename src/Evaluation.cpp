/* ------------------------------------------------------------------------
*
 * SPDX-License-Identifier: GPL-3-or-later
 * Copyright (C) 2023 - 2024 by Matteo Malvestiti
 *
 * This file is part of ion_propulsion.
 *
 * ion_propulsion is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ion_propulsion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ion_propulsion; see the file COPYING.  If not, see
 * <https://www.gnu.org/licenses/>.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Matteo Malvestiti, Politecnico di Milano, 2024
 *
 */

#include "Evaluation.h"

#include <Constants.h>

namespace IonPropulsion{
	using namespace dealii;
	namespace Evaluation{
		// ------------------------------------------------------
		// EvaluationBase
		// ------------------------------------------------------
		template <int dim>
		EvaluationBase<dim>::EvaluationBase(const unsigned degree):
		mapping(degree)
		, mpi_communicator(MPI_COMM_WORLD)
		, n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
		, this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
		, pcout(std::cout, (this_mpi_process == 0)){};

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
		const unsigned degree,
		const Point<dim> &evaluation_point)
			: EvaluationBase<dim>(degree), evaluation_point(evaluation_point)
		{}

		template <int dim>
		std::pair<std::string, double>
		PointValueEvaluation<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const PETScWrappers::MPI::Vector & locally_relevant_solution,
			const parallel::distributed::Triangulation<dim> &       triangulation) const
		{
			(void)triangulation;
			double point_value = 1e20;
			double min_distance = std::numeric_limits<double>::max();
			double current_distance = std::numeric_limits<double>::max();
			typename DoFHandler<dim>::active_cell_iterator nearest_cell;
			unsigned int nearest_vertex = 0;

			for (const auto &cell : dof_handler.active_cell_iterators()) {
				if (cell->is_locally_owned()){
					for (const auto vertex : cell->vertex_indices()) {
						current_distance = cell->vertex(vertex).distance(evaluation_point);
						if (current_distance < min_distance){
							min_distance = current_distance;
							nearest_cell = cell;
							nearest_vertex = vertex;
						}
					}
				}
			}

			// Compute the global minimum distance across all MPI ranks
			double overall_minimum_distance=std::numeric_limits<double>::max();
			MPI_Allreduce(&min_distance, &overall_minimum_distance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

			// 1. Identify the sending rank (rank with the closest point)
			int sending_rank = -1;
			if (std::fabs(min_distance - overall_minimum_distance) < 1.e-16) {
				sending_rank = this->this_mpi_process;
			}

			// 2. Communicate the sending rank ID to all ranks
			int global_sending_rank;
			MPI_Allreduce(&sending_rank, &global_sending_rank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			// 3. If this rank is the sending rank, send the data to rank 0
			if (this->this_mpi_process == global_sending_rank) {
				point_value = locally_relevant_solution(nearest_cell->vertex_dof_index(nearest_vertex, 0));
				cout<<"   Closest point found at distance "<<min_distance
						<<" from the required point by rank: "<<this->this_mpi_process<<std::endl;

				double exact_error = std::fabs(point_value-EXACT_POINT_VALUE);

				std::pair<std::string, double> pair = std::make_pair("ex POINT err",exact_error);

				if (this->this_mpi_process == 0) {
					return pair;
				} else {
					int key_length = pair.first.length();
					MPI_Send(&key_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // Send string length first
					MPI_Send(pair.first.c_str(), key_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD); // Send string
					MPI_Send(&pair.second, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // Send double value
				}
			}

			// 4. If this rank is rank 0, receive data from the sending rank
			if (this->this_mpi_process == 0) {
				if (global_sending_rank != 0) {
					int key_length = 0;
					MPI_Recv(&key_length, 1, MPI_INT, global_sending_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive string length
					std::vector<char> key_buffer(key_length + 1, '\0');
					MPI_Recv(key_buffer.data(), key_length, MPI_CHAR, global_sending_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive string
					double value;
					MPI_Recv(&value, 1, MPI_DOUBLE, global_sending_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive double value
					std::string key(key_buffer.data());

					return std::make_pair(key, value);
				}
			}

			// If this rank neither finds the point nor is rank 0, do nothing
			return std::make_pair("null", -1.0e-20);
		}

		// ------------------------------------------------------
		// FluxEvaluation
		// ------------------------------------------------------

		template <int dim>
		FluxEvaluation<dim>::FluxEvaluation(const unsigned degree, const std::set<unsigned int> &boundary_ids):
		EvaluationBase<dim>(degree), boundary_ids(boundary_ids)
		{}

		template <int dim>
		std::pair<std::string, double>
		FluxEvaluation<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const PETScWrappers::MPI::Vector & solution,
			const parallel::distributed::Triangulation<dim> &       triangulation) const
		{
			(void)triangulation;
			double flux = 0;
			const QGauss<dim-1> face_quadrature(dof_handler.get_fe().degree + 1);
			FEFaceValues<dim> fe_face_values(this->mapping, dof_handler.get_fe(),
																			 face_quadrature,
																			 update_gradients | update_normal_vectors | update_JxW_values);
			const unsigned int n_face_q_points = face_quadrature.size();
			std::vector<Tensor<1, dim>> solution_gradients(n_face_q_points);

			for (const auto &cell : dof_handler.active_cell_iterators())
				if (cell->is_locally_owned()){
					for (const auto &face : cell->face_iterators())
						if (face->at_boundary() && boundary_ids.count(face->boundary_id()) > 0 ) {
							fe_face_values.reinit(cell, face);
							fe_face_values.get_function_gradients(solution, solution_gradients);
							for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) {
								const Tensor<1, dim> &n = fe_face_values.normal_vector(q_point);
								flux += ((-solution_gradients[q_point]) * (-n)) * fe_face_values.JxW(q_point);
							}
						}
				}

			MPI_Barrier(this->mpi_communicator);
			double global_flux = 0.0;
			MPI_Reduce(&flux, &global_flux, 1, MPI_DOUBLE, MPI_SUM, 0, this->mpi_communicator);

			if (LOAD_FROM_SETUP != 0 && LOAD_FROM_SETUP != 11) {
				double exact_error = std::fabs(global_flux-EXACT_FLUX);
				return std::make_pair("std FLUX err",exact_error);
			} else {
				return std::make_pair("std FLUX value",global_flux);
			}
		}

		// ------------------------------------------------------
		// GridOutput
		// ------------------------------------------------------

		template <int dim>
		GridOutput<dim>::GridOutput(const unsigned degree, const std::string &output_name_base)
			: EvaluationBase<dim>(degree), output_name_base(output_name_base)
		{}


		template <int dim>
		std::pair<std::string, double>
		GridOutput<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const PETScWrappers::MPI::Vector & /*solution*/,
			const parallel::distributed::Triangulation<dim> &       triangulation) const
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
		L2_error_estimate<dim>::L2_error_estimate(const unsigned degree, const Function<dim> & analytical_solution)
			: EvaluationBase<dim>(degree), analytical_solution(&analytical_solution)
		{}


		template <int dim>
		std::pair<std::string, double>
		L2_error_estimate<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const PETScWrappers::MPI::Vector & locally_relevant_solution,
			const parallel::distributed::Triangulation<dim> &       triangulation) const
		{
			// Set up an error vector
			Vector<double> difference_per_cell(triangulation.n_active_cells());

			// Compute the difference between the finite element solution and the exact solution
			const QGauss<dim> quadrature_formula(2*dof_handler.get_fe().degree + 1);
			VectorTools::integrate_difference(
					this->mapping,
					dof_handler,
					locally_relevant_solution,
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
		H1_error_estimate<dim>::H1_error_estimate(const unsigned degree, const Function<dim> & analytical_solution)
			: EvaluationBase<dim>(degree), analytical_solution(&analytical_solution)
		{}


		template <int dim>
		std::pair<std::string, double>
		H1_error_estimate<dim>::operator()(
			const DoFHandler<dim> &dof_handler,
			const PETScWrappers::MPI::Vector & locally_relevant_solution,
			const parallel::distributed::Triangulation<dim> &       triangulation) const
		{
			// Set up an error vector
			Vector<double> difference_per_cell(triangulation.n_active_cells());

			// Compute the difference between the finite element solution and the exact solution
			const QGauss<dim> quadrature_formula(2*dof_handler.get_fe().degree + 1);
			VectorTools::integrate_difference(
					this->mapping,
					dof_handler,
					locally_relevant_solution,
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