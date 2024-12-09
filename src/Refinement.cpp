#include "Refinement.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{
    // ------------------------------------------------------
    // RefinementGlobal
    // ------------------------------------------------------
    template <int dim>
    RefinementGlobal<dim>::RefinementGlobal(
      Triangulation<dim> &       coarse_grid,
      const FiniteElement<dim> & fe,
      const Quadrature<dim> &    quadrature,
      const Quadrature<dim - 1> &face_quadrature,
      const Function<dim> &      rhs_function,
      const Function<dim> &      boundary_values)
      : Base<dim>(coarse_grid)
      , PrimalSolver<dim>(coarse_grid,
                          fe,
                          quadrature,
                          face_quadrature,
                          rhs_function,
                          boundary_values)
      {}

    template <int dim>
    void RefinementGlobal<dim>::refine_grid()
    {
      Vector<double> old_Rg_values = this->Rg_vector;

      SolutionTransfer<dim> solution_transfer(this->dof_handler);
      solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_values);

      this->triangulation->refine_global(1);

      this->dof_handler.distribute_dofs(*(this->fe));

      this->Rg_vector.reinit(this->dof_handler.n_dofs());
      solution_transfer.interpolate(old_Rg_values, this->Rg_vector);

      this->construct_Rg_vector();
    }

    // ------------------------------------------------------
    // RefinementKelly
    // ------------------------------------------------------
    template <int dim>
    RefinementKelly<dim>::RefinementKelly(
      Triangulation<dim> &       coarse_grid,
      const FiniteElement<dim> & fe,
      const Quadrature<dim> &    quadrature,
      const Quadrature<dim - 1> &face_quadrature,
      const Function<dim> &      rhs_function,
      const Function<dim> &      boundary_values)
      : Base<dim>(coarse_grid)
      , PrimalSolver<dim>(coarse_grid,
                          fe,
                          quadrature,
                          face_quadrature,
                          rhs_function,
                          boundary_values)
    {}



    template <int dim>
    void RefinementKelly<dim>::refine_grid()
    {
      Vector<float> estimated_error_per_cell(
        this->triangulation->n_active_cells());
      KellyErrorEstimator<dim>::estimate(
        this->dof_handler,
        QGauss<dim - 1>(this->fe->degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        this->solution,
        estimated_error_per_cell);
      GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
                                                      estimated_error_per_cell,
                                                      0.3,
                                                      0.03);
      this->triangulation->execute_coarsening_and_refinement();
    }

    // ------------------------------------------------------
    // RefinementWeightedKelly
    // ------------------------------------------------------
    template <int dim>
    RefinementWeightedKelly<dim>::RefinementWeightedKelly(
      Triangulation<dim> &       coarse_grid,
      const FiniteElement<dim> & fe,
      const Quadrature<dim> &    quadrature,
      const Quadrature<dim - 1> &face_quadrature,
      const Function<dim> &      rhs_function,
      const Function<dim> &      boundary_values,
      const Function<dim> &      weighting_function)
      : Base<dim>(coarse_grid)
      , PrimalSolver<dim>(coarse_grid,
                          fe,
                          quadrature,
                          face_quadrature,
                          rhs_function,
                          boundary_values)
      , weighting_function(&weighting_function)
    {}

    template <int dim>
    void RefinementWeightedKelly<dim>::refine_grid()
    {
      Vector<float> estimated_error_per_cell(
        this->triangulation->n_active_cells());
      std::map<types::boundary_id, const Function<dim> *> dummy_function_map;
      KellyErrorEstimator<dim>::estimate(this->dof_handler,
                                         *this->face_quadrature,
                                         dummy_function_map,
                                         this->solution,
                                         estimated_error_per_cell);

      for (const auto &cell : this->dof_handler.active_cell_iterators())
        estimated_error_per_cell(cell->active_cell_index()) *=
          weighting_function->value(cell->center());

      GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
                                                      estimated_error_per_cell,
                                                      0.3,
                                                      0.03);
      this->triangulation->execute_coarsening_and_refinement();
    }

    template <int dim>
    WeightedResidual<dim>::CellData::CellData(
      const FiniteElement<dim> &fe,
      const Quadrature<dim> &   quadrature,
      const Function<dim> &     right_hand_side)
      : fe_values(fe,
                  quadrature,
                  update_values | update_hessians | update_quadrature_points |
                    update_JxW_values)
      , right_hand_side(&right_hand_side)
      , cell_residual(quadrature.size())
      , rhs_values(quadrature.size())
      , dual_weights(quadrature.size())
      , cell_laplacians(quadrature.size())
    {}

    template <int dim>
    WeightedResidual<dim>::CellData::CellData(const CellData &cell_data)
      : fe_values(cell_data.fe_values.get_fe(),
                  cell_data.fe_values.get_quadrature(),
                  update_values | update_hessians | update_quadrature_points |
                    update_JxW_values)
      , right_hand_side(cell_data.right_hand_side)
      , cell_residual(cell_data.cell_residual)
      , rhs_values(cell_data.rhs_values)
      , dual_weights(cell_data.dual_weights)
      , cell_laplacians(cell_data.cell_laplacians)
    {}

    template <int dim>
    WeightedResidual<dim>::FaceData::FaceData(
      const FiniteElement<dim> & fe,
      const Quadrature<dim - 1> &face_quadrature)
      : fe_face_values_cell(fe,
                            face_quadrature,
                            update_values | update_gradients |
                              update_JxW_values | update_normal_vectors)
      , fe_face_values_neighbor(fe,
                                face_quadrature,
                                update_values | update_gradients |
                                  update_JxW_values | update_normal_vectors)
      , fe_subface_values_cell(fe, face_quadrature, update_gradients)
    {
      const unsigned int n_face_q_points = face_quadrature.size();

      jump_residual.resize(n_face_q_points);
      dual_weights.resize(n_face_q_points);
      cell_grads.resize(n_face_q_points);
      neighbor_grads.resize(n_face_q_points);
    }

    template <int dim>
    WeightedResidual<dim>::FaceData::FaceData(const FaceData &face_data)
      : fe_face_values_cell(face_data.fe_face_values_cell.get_fe(),
                            face_data.fe_face_values_cell.get_quadrature(),
                            update_values | update_gradients |
                              update_JxW_values | update_normal_vectors)
      , fe_face_values_neighbor(
          face_data.fe_face_values_neighbor.get_fe(),
          face_data.fe_face_values_neighbor.get_quadrature(),
          update_values | update_gradients | update_JxW_values |
            update_normal_vectors)
      , fe_subface_values_cell(
          face_data.fe_subface_values_cell.get_fe(),
          face_data.fe_subface_values_cell.get_quadrature(),
          update_gradients)
      , jump_residual(face_data.jump_residual)
      , dual_weights(face_data.dual_weights)
      , cell_grads(face_data.cell_grads)
      , neighbor_grads(face_data.neighbor_grads)
    {}

    template <int dim>
    WeightedResidual<dim>::WeightedResidualScratchData::
      WeightedResidualScratchData(
        const FiniteElement<dim> & primal_fe,
        const Quadrature<dim> &    primal_quadrature,
        const Quadrature<dim - 1> &primal_face_quadrature,
        const Function<dim> &      rhs_function,
        const Vector<double> &     primal_solution,
        const Vector<double> &     dual_weights)
      : cell_data(primal_fe, primal_quadrature, rhs_function)
      , face_data(primal_fe, primal_face_quadrature)
      , primal_solution(primal_solution)
      , dual_weights(dual_weights)
    {}

    template <int dim>
    WeightedResidual<dim>::WeightedResidualScratchData::
      WeightedResidualScratchData(
        const WeightedResidualScratchData &scratch_data)
      : cell_data(scratch_data.cell_data)
      , face_data(scratch_data.face_data)
      , primal_solution(scratch_data.primal_solution)
      , dual_weights(scratch_data.dual_weights)
    {}

    template <int dim>
    WeightedResidual<dim>::WeightedResidual(
      Triangulation<dim> &                           coarse_grid,
      const FiniteElement<dim> &                     primal_fe,
      const FiniteElement<dim> &                     dual_fe,
      const Quadrature<dim> &                        quadrature,
      const Quadrature<dim - 1> &                    face_quadrature,
      const Function<dim> &                          rhs_function,
      const Function<dim> &                          bv,
      const DualFunctional::DualFunctionalBase<dim> &dual_functional)
      : Base<dim>(coarse_grid)
      , PrimalSolver<dim>(coarse_grid,
                          primal_fe,
                          quadrature,
                          face_quadrature,
                          rhs_function,
                          bv)
      , DualSolver<dim>(coarse_grid,
                        dual_fe,
                        quadrature,
                        face_quadrature,
                        dual_functional)
    {}

    template <int dim>
    void WeightedResidual<dim>::solve_problem()
    {
      Threads::TaskGroup<void> tasks;
      tasks +=
        Threads::new_task(&WeightedResidual<dim>::solve_primal_problem, *this);
      tasks +=
        Threads::new_task(&WeightedResidual<dim>::solve_dual_problem, *this);
      tasks.join_all();
    }


    template <int dim>
    void WeightedResidual<dim>::solve_primal_problem()
    {
      PrimalSolver<dim>::solve_problem();
    }

    template <int dim>
    void WeightedResidual<dim>::solve_dual_problem()
    {
      DualSolver<dim>::solve_problem();
    }


    template <int dim>
    void WeightedResidual<dim>::postprocess(
      const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      PrimalSolver<dim>::postprocess(postprocessor);
    }


    template <int dim>
    unsigned int WeightedResidual<dim>::n_dofs() const
    {
      return PrimalSolver<dim>::n_dofs();
    }

    template <int dim>
    void WeightedResidual<dim>::refine_grid()
    {
      Vector<float> error_indicators(this->triangulation->n_active_cells());
      estimate_error(error_indicators);

      for (float &error_indicator : error_indicators)
        error_indicator = std::fabs(error_indicator);

      GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                        error_indicators,
                                                        0.8,
                                                        0.02);
      this->triangulation->execute_coarsening_and_refinement();
    }

    template <int dim>
    void WeightedResidual<dim>::output_solution() const
    {
      AffineConstraints<double> primal_hanging_node_constraints;
      DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler,
                                              primal_hanging_node_constraints);
      primal_hanging_node_constraints.close();
      Vector<double> dual_solution(PrimalSolver<dim>::dof_handler.n_dofs());
      FETools::interpolate(DualSolver<dim>::dof_handler,
                           DualSolver<dim>::solution,
                           PrimalSolver<dim>::dof_handler,
                           primal_hanging_node_constraints,
                           dual_solution);

      DataOut<dim> data_out;
      data_out.attach_dof_handler(PrimalSolver<dim>::dof_handler);

      // Add the data vectors for which we want output. Add them both, the
      // <code>DataOut</code> functions can handle as many data vectors as you
      // wish to write to output:
      data_out.add_data_vector(PrimalSolver<dim>::solution, "primal_solution");
      data_out.add_data_vector(dual_solution, "dual_solution");

      data_out.build_patches();

      std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                        ".vtu");
      data_out.write(out, DataOutBase::vtu);
    }

    template <int dim>
    void
    WeightedResidual<dim>::estimate_error(Vector<float> &error_indicators) const
    {

      AffineConstraints<double> dual_hanging_node_constraints;
      DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler,
                                              dual_hanging_node_constraints);
      dual_hanging_node_constraints.close();
      Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs());
      FETools::interpolate(PrimalSolver<dim>::dof_handler,
                           PrimalSolver<dim>::solution,
                           DualSolver<dim>::dof_handler,
                           dual_hanging_node_constraints,
                           primal_solution);


      AffineConstraints<double> primal_hanging_node_constraints;
      DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler,
                                              primal_hanging_node_constraints);
      primal_hanging_node_constraints.close();
      Vector<double> dual_weights(DualSolver<dim>::dof_handler.n_dofs());
      FETools::interpolation_difference(DualSolver<dim>::dof_handler,
                                        dual_hanging_node_constraints,
                                        DualSolver<dim>::solution,
                                        PrimalSolver<dim>::dof_handler,
                                        primal_hanging_node_constraints,
                                        dual_weights);

      FaceIntegrals face_integrals;
      for (const auto &cell :
           DualSolver<dim>::dof_handler.active_cell_iterators())
        for (const auto &face : cell->face_iterators())
          face_integrals[face] = -1e20;

      auto worker = [this,
                     &error_indicators,
                     &face_integrals](const active_cell_iterator & cell,
                                      WeightedResidualScratchData &scratch_data,
                                      WeightedResidualCopyData &   copy_data) {
        this->estimate_on_one_cell(
          cell, scratch_data, copy_data, error_indicators, face_integrals);
      };

      auto do_nothing_copier =
        std::function<void(const WeightedResidualCopyData &)>();

      // Then hand it all off to WorkStream::run() to compute the
      // estimators for all cells in parallel:
      WorkStream::run(
        DualSolver<dim>::dof_handler.begin_active(),
        DualSolver<dim>::dof_handler.end(),
        worker,
        do_nothing_copier,
        WeightedResidualScratchData(*DualSolver<dim>::fe,
                                    *DualSolver<dim>::quadrature,
                                    *DualSolver<dim>::face_quadrature,
                                    *this->rhs_function,
                                    primal_solution,
                                    dual_weights),
        WeightedResidualCopyData());

      unsigned int present_cell = 0;
      for (const auto &cell :
           DualSolver<dim>::dof_handler.active_cell_iterators())
        {
          for (const auto &face : cell->face_iterators())
            {
              Assert(face_integrals.find(face) != face_integrals.end(),
                     ExcInternalError());
              error_indicators(present_cell) -= 0.5 * face_integrals[face];
            }
          ++present_cell;
        }
      std::cout << "   Estimated error="
                << std::accumulate(error_indicators.begin(),
                                   error_indicators.end(),
                                   0.)
                << std::endl;
    }

    template <int dim>
    void WeightedResidual<dim>::estimate_on_one_cell(
      const active_cell_iterator & cell,
      WeightedResidualScratchData &scratch_data,
      WeightedResidualCopyData &   copy_data,
      Vector<float> &              error_indicators,
      FaceIntegrals &              face_integrals) const
    {
      (void)copy_data;
      integrate_over_cell(cell,
                          scratch_data.primal_solution,
                          scratch_data.dual_weights,
                          scratch_data.cell_data,
                          error_indicators);
      for (const auto face_no : cell->face_indices())
        {
          if (cell->face(face_no)->at_boundary())
            {
              face_integrals[cell->face(face_no)] = 0;
              continue;
            }
          if ((cell->neighbor(face_no)->has_children() == false) &&
              (cell->neighbor(face_no)->level() == cell->level()) &&
              (cell->neighbor(face_no)->index() < cell->index()))
            continue;
          if (cell->at_boundary(face_no) == false)
            if (cell->neighbor(face_no)->level() < cell->level())
              continue;
          if (cell->face(face_no)->has_children() == false)
            integrate_over_regular_face(cell,
                                        face_no,
                                        scratch_data.primal_solution,
                                        scratch_data.dual_weights,
                                        scratch_data.face_data,
                                        face_integrals);
          else
            integrate_over_irregular_face(cell,
                                          face_no,
                                          scratch_data.primal_solution,
                                          scratch_data.dual_weights,
                                          scratch_data.face_data,
                                          face_integrals);
        }
    }

    template <int dim>
    void WeightedResidual<dim>::integrate_over_cell(
      const active_cell_iterator &cell,
      const Vector<double> &      primal_solution,
      const Vector<double> &      dual_weights,
      CellData &                  cell_data,
      Vector<float> &             error_indicators) const
    {
      // The tasks to be done are what appears natural from looking at the
      // error estimation formula: first get the right hand side and Laplacian
      // of the numerical solution at the quadrature points for the cell
      // residual,
      cell_data.fe_values.reinit(cell);
      cell_data.right_hand_side->value_list(
        cell_data.fe_values.get_quadrature_points(), cell_data.rhs_values);
      cell_data.fe_values.get_function_laplacians(primal_solution,
                                                  cell_data.cell_laplacians);

      // ...then get the dual weights...
      cell_data.fe_values.get_function_values(dual_weights,
                                              cell_data.dual_weights);

      // ...and finally build the sum over all quadrature points and store it
      // with the present cell:
      double sum = 0;
      for (unsigned int p = 0; p < cell_data.fe_values.n_quadrature_points; ++p)
        sum += ((cell_data.rhs_values[p] + cell_data.cell_laplacians[p]) *
                cell_data.dual_weights[p] * cell_data.fe_values.JxW(p));
      error_indicators(cell->active_cell_index()) += sum;
    }


    // @sect4{Computing edge term error contributions -- 1}

    // On the other hand, computation of the edge terms for the error estimate
    // is not so simple. First, we have to distinguish between faces with and
    // without hanging nodes. Because it is the simple case, we first consider
    // the case without hanging nodes on a face (let's call this the `regular'
    // case):
    template <int dim>
    void WeightedResidual<dim>::integrate_over_regular_face(
      const active_cell_iterator &cell,
      const unsigned int          face_no,
      const Vector<double> &      primal_solution,
      const Vector<double> &      dual_weights,
      FaceData &                  face_data,
      FaceIntegrals &             face_integrals) const
    {
      const unsigned int n_q_points =
        face_data.fe_face_values_cell.n_quadrature_points;

      // The first step is to get the values of the gradients at the
      // quadrature points of the finite element field on the present
      // cell. For this, initialize the <code>FEFaceValues</code> object
      // corresponding to this side of the face, and extract the gradients
      // using that object.
      face_data.fe_face_values_cell.reinit(cell, face_no);
      face_data.fe_face_values_cell.get_function_gradients(
        primal_solution, face_data.cell_grads);

      // The second step is then to extract the gradients of the finite
      // element solution at the quadrature points on the other side of the
      // face, i.e. from the neighboring cell.
      //
      // For this, do a sanity check before: make sure that the neighbor
      // actually exists (yes, we should not have come here if the neighbor
      // did not exist, but in complicated software there are bugs, so better
      // check this), and if this is not the case throw an error.
      Assert(cell->neighbor(face_no).state() == IteratorState::valid,
             ExcInternalError());
      // If we have that, then we need to find out with which face of the
      // neighboring cell we have to work, i.e. the <code>how-many'th</code> the
      // neighbor the present cell is of the cell behind the present face. For
      // this, there is a function, and we put the result into a variable with
      // the name <code>neighbor_neighbor</code>:
      const unsigned int neighbor_neighbor =
        cell->neighbor_of_neighbor(face_no);
      // Then define an abbreviation for the neighbor cell, initialize the
      // <code>FEFaceValues</code> object on that cell, and extract the
      // gradients on that cell:
      const active_cell_iterator neighbor = cell->neighbor(face_no);
      face_data.fe_face_values_neighbor.reinit(neighbor, neighbor_neighbor);
      face_data.fe_face_values_neighbor.get_function_gradients(
        primal_solution, face_data.neighbor_grads);

      // Now that we have the gradients on this and the neighboring cell,
      // compute the jump residual by multiplying the jump in the gradient
      // with the normal vector:
      for (unsigned int p = 0; p < n_q_points; ++p)
        face_data.jump_residual[p] =
          ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) *
           face_data.fe_face_values_cell.normal_vector(p));

      // Next get the dual weights for this face:
      face_data.fe_face_values_cell.get_function_values(dual_weights,
                                                        face_data.dual_weights);

      // Finally, we have to compute the sum over jump residuals, dual
      // weights, and quadrature weights, to get the result for this face:
      double face_integral = 0;
      for (unsigned int p = 0; p < n_q_points; ++p)
        face_integral +=
          (face_data.jump_residual[p] * face_data.dual_weights[p] *
           face_data.fe_face_values_cell.JxW(p));


      Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(),
             ExcInternalError());
      Assert(face_integrals[cell->face(face_no)] == -1e20, ExcInternalError());

      face_integrals[cell->face(face_no)] = face_integral;
    }

    template <int dim>
    void WeightedResidual<dim>::integrate_over_irregular_face(
      const active_cell_iterator &cell,
      const unsigned int          face_no,
      const Vector<double> &      primal_solution,
      const Vector<double> &      dual_weights,
      FaceData &                  face_data,
      FaceIntegrals &             face_integrals) const
    {
      const unsigned int n_q_points =
        face_data.fe_face_values_cell.n_quadrature_points;

      const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
      const typename DoFHandler<dim>::cell_iterator neighbor =
        cell->neighbor(face_no);
      Assert(neighbor.state() == IteratorState::valid, ExcInternalError());
      Assert(neighbor->has_children(), ExcInternalError());
      (void)neighbor;

      const unsigned int neighbor_neighbor =
        cell->neighbor_of_neighbor(face_no);

      // Then simply do everything we did in the previous function for one
      // face for all the sub-faces now:
      for (unsigned int subface_no = 0; subface_no < face->n_children();
           ++subface_no)
        {

          const active_cell_iterator neighbor_child =
            cell->neighbor_child_on_subface(face_no, subface_no);
          Assert(neighbor_child->face(neighbor_neighbor) ==
                   cell->face(face_no)->child(subface_no),
                 ExcInternalError());

          // Now start the work by again getting the gradient of the solution
          // first at this side of the interface,
          face_data.fe_subface_values_cell.reinit(cell, face_no, subface_no);
          face_data.fe_subface_values_cell.get_function_gradients(
            primal_solution, face_data.cell_grads);
          // then at the other side,
          face_data.fe_face_values_neighbor.reinit(neighbor_child,
                                                   neighbor_neighbor);
          face_data.fe_face_values_neighbor.get_function_gradients(
            primal_solution, face_data.neighbor_grads);

          for (unsigned int p = 0; p < n_q_points; ++p)
            face_data.jump_residual[p] =
              ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) *
               face_data.fe_face_values_neighbor.normal_vector(p));

          // Then get dual weights:
          face_data.fe_face_values_neighbor.get_function_values(
            dual_weights, face_data.dual_weights);

          // At last, sum up the contribution of this sub-face, and set it in
          // the global map:
          double face_integral = 0;
          for (unsigned int p = 0; p < n_q_points; ++p)
            face_integral +=
              (face_data.jump_residual[p] * face_data.dual_weights[p] *
               face_data.fe_face_values_neighbor.JxW(p));
          face_integrals[neighbor_child->face(neighbor_neighbor)] =
            face_integral;
        }

      double sum = 0;
      for (unsigned int subface_no = 0; subface_no < face->n_children();
           ++subface_no)
        {
          Assert(face_integrals.find(face->child(subface_no)) !=
                   face_integrals.end(),
                 ExcInternalError());
          Assert(face_integrals[face->child(subface_no)] != -1e20,
                 ExcInternalError());

          sum += face_integrals[face->child(subface_no)];
        }
      // Finally store the value with the parent face.
      face_integrals[face] = sum;
    }
    // Template instantiation
    template class RefinementGlobal<2>;
    template class RefinementKelly<2>;
    template class RefinementWeightedKelly<2>;
    template class WeightedResidual<2>;
  } // namespace LaplaceSolver
} // namespace IonPropulsion