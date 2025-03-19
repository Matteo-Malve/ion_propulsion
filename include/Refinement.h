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

#ifndef REFINEMENT_H
#define REFINEMENT_H
#include "includes.h"
#include "LaplaceSolver.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{
    // ------------------------------------------------------
    // RefinementGlobal
    // ------------------------------------------------------
    template <int dim>
    class RefinementGlobal : public PrimalSolver<dim>
    {
    public:
      RefinementGlobal(parallel::distributed::Triangulation<dim> &       coarse_grid,
                       const FiniteElement<dim> & fe,
                       const Quadrature<dim> &    quadrature,
                       const Quadrature<dim - 1> &face_quadrature,
                       const Function<dim> &      rhs_function,
                       const Function<dim> &      boundary_values,
                       const unsigned degree);

      virtual void refine_grid() override;
      void print_convergence_table() const override;
    };

    // ------------------------------------------------------
    // WeightedResidual
    // ------------------------------------------------------

    template <int dim>
    class WeightedResidual : public PrimalSolver<dim>, public DualSolver<dim>
    {
    public:
      WeightedResidual(
        parallel::distributed::Triangulation<dim> &                           coarse_grid,
        const FiniteElement<dim> &                     primal_fe,
        const FiniteElement<dim> &                     dual_fe,
        const Quadrature<dim> &                        quadrature,
        const Quadrature<dim - 1> &                    face_quadrature,
        const Function<dim> &                          rhs_function,
        const Function<dim> &                          boundary_values,
        const DualFunctional::DualFunctionalBase<dim> &dual_functional,
        const unsigned degree);

      virtual void solve_problem() override;

      virtual void postprocess(
        const Evaluation::EvaluationBase<dim> &postprocessor) const override;

      virtual unsigned int n_dofs() const override;

      virtual void refine_grid() override;

      virtual void output_solution()  override;

      void update_convergence_table() override;
      void print_convergence_table() const override;

    private:

      void solve_primal_problem();
      void solve_dual_problem();

      using active_cell_iterator =
        typename DoFHandler<dim>::active_cell_iterator;

      using FaceIntegrals =
        typename std::map<typename DoFHandler<dim>::face_iterator, double>;

      struct CellData
      {
        MappingQ<dim>  mapping;
        FEValues<dim>                           fe_values;
        const SmartPointer<const Function<dim>> right_hand_side;

        std::vector<double> cell_residual;
        std::vector<double> rhs_values;
        std::vector<double> dual_weights;
        std::vector<double> cell_laplacians;
        CellData(const FiniteElement<dim> &fe,
                 const Quadrature<dim> &   quadrature,
                 const Function<dim> &     right_hand_side,
                 const MappingQ<dim> & mapping);
        CellData(const CellData &cell_data);
      };

      struct FaceData
      {
        MappingQ<dim>  mapping;
        FEFaceValues<dim>    fe_face_values_cell;
        FEFaceValues<dim>    fe_face_values_neighbor;
        FESubfaceValues<dim> fe_subface_values_cell;
        FESubfaceValues<dim> fe_subface_values_neighbor;

        std::vector<double>                  jump_residual;
        std::vector<double>                  dual_weights;
        typename std::vector<Tensor<1, dim>> cell_grads;
        typename std::vector<Tensor<1, dim>> neighbor_grads;
        FaceData(const FiniteElement<dim> & fe,
                 const Quadrature<dim - 1> &face_quadrature,
                 const MappingQ<dim> & mapping);
        FaceData(const FaceData &face_data);
      };

      struct WeightedResidualScratchData
      {

        WeightedResidualScratchData(
          const FiniteElement<dim> & primal_fe,
          const Quadrature<dim> &    primal_quadrature,
          const Quadrature<dim - 1> &primal_face_quadrature,
          const Function<dim> &      rhs_function,
          const PETScWrappers::MPI::Vector &     primal_solution,
          const PETScWrappers::MPI::Vector &     dual_weights,
          const MappingQ<dim> & mapping);

        WeightedResidualScratchData(
          const WeightedResidualScratchData &scratch_data);

        MappingQ<dim>  mapping;
        CellData       cell_data;
        FaceData       face_data;
        PETScWrappers::MPI::Vector primal_solution;
        PETScWrappers::MPI::Vector dual_weights;
      };

      struct WeightedResidualCopyData
      {};

      void estimate_error(Vector<float> &error_indicators) const;

      void estimate_on_one_cell(const active_cell_iterator & cell,
                                WeightedResidualScratchData &scratch_data,
                                WeightedResidualCopyData &   copy_data,
                                Vector<float> &              error_indicators,
                                FaceIntegrals &face_integrals) const;

      void integrate_over_cell(const active_cell_iterator &cell,
                               const PETScWrappers::MPI::Vector &      primal_solution,
                               const PETScWrappers::MPI::Vector &      dual_weights,
                               CellData &                  cell_data,
                               Vector<float> &error_indicators) const;

      void integrate_over_regular_face(const active_cell_iterator &cell,
                                       const unsigned int          face_no,
                                       const PETScWrappers::MPI::Vector &primal_solution,
                                       const PETScWrappers::MPI::Vector &dual_weights,
                                       FaceData &            face_data,
                                       FaceIntegrals &face_integrals) const;
      void integrate_over_irregular_face(const active_cell_iterator &cell,
                                         const unsigned int          face_no,
                                         const PETScWrappers::MPI::Vector &primal_solution,
                                         const PETScWrappers::MPI::Vector &dual_weights,
                                         FaceData &            face_data,
                                         FaceIntegrals &face_integrals) const;
      void integrate_over_irregular_face_flipped(const active_cell_iterator &cell,
                                         const unsigned int          face_no,
                                         const PETScWrappers::MPI::Vector &primal_solution,
                                         const PETScWrappers::MPI::Vector &dual_weights,
                                         FaceData &            face_data,
                                         FaceIntegrals &face_integrals) const;
    };



  } // namespace LaplaceSolver
} // namespace IonPropulsion




#endif