#include "Refinement.h"

#include <Constants.h>

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{
    // ------------------------------------------------------
    // RefinementGlobal
    // ------------------------------------------------------
    template <int dim>
    RefinementGlobal<dim>::RefinementGlobal(
      parallel::distributed::Triangulation<dim> &       coarse_grid,
      const FiniteElement<dim> & fe,
      const Quadrature<dim> &    quadrature,
      const Quadrature<dim - 1> &face_quadrature,
      const Function<dim> &      rhs_function,
      const Function<dim> &      boundary_values,
      const unsigned degree)
      : Base<dim>(coarse_grid)
      , PrimalSolver<dim>(coarse_grid,
                          fe,
                          quadrature,
                          face_quadrature,
                          rhs_function,
                          boundary_values,
                          degree)
      {}

    template <int dim>
    void RefinementGlobal<dim>::refine_grid()
    {
      TimerOutput::Scope t(this->computing_timer, "refine");
      // REFINE EVERYWHERE
      if (REFINEMENT_CRITERION == 1) {
        if (MANUAL_LIFTING_ON) {    // TODO
          /*Vector<double> old_Rg_values = this->locally_relevantRg_vector;
          SolutionTransfer<dim> solution_transfer(this->dof_handler);
          solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_values);
          this->triangulation->refine_global(1);
          this->dof_handler.distribute_dofs(*(this->fe));
          this->locally_relevantRg_vector.reinit(this->dof_handler.n_dofs());
          solution_transfer.interpolate(old_Rg_values, this->locally_relevantRg_vector);
          this->construct_Rg_vector();*/
        } else {
          this->triangulation->refine_global(1);
        }
      }

      // ONLY REFINE AROUND EMITTER
      else if (REFINEMENT_CRITERION == 4) {

        std::unique_ptr<const std::set<unsigned int>> emitter_boundary_ids_set_ptr;

        if (LOAD_FROM_SETUP==0)
          emitter_boundary_ids_set_ptr = std::make_unique<const std::set<unsigned int>>(std::set<unsigned int>{1});
        else  if (LOAD_FROM_SETUP==11)
          emitter_boundary_ids_set_ptr = std::make_unique<const std::set<unsigned int>>(std::set<unsigned int>{1,2});
        else
          emitter_boundary_ids_set_ptr = std::make_unique<const std::set<unsigned int>>(std::set<unsigned int>{1});

        Vector<float> criteria(this->triangulation->n_active_cells());
        for (auto &cell : this->triangulation->active_cell_iterators()) {
          if (cell->is_locally_owned())
            if (cell->at_boundary())
              for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
                if (cell->face(face)->at_boundary() && emitter_boundary_ids_set_ptr->count(cell->face(face)->boundary_id()) > 0) {
                  cell->set_refine_flag();
                  break;
                }
        }

        if (MANUAL_LIFTING_ON) {    //TODO
          /*this->triangulation->prepare_coarsening_and_refinement();
          SolutionTransfer<dim> solution_transfer(this->dof_handler);

          Vector<double> old_Rg_values = this->locally_relevantRg_vector;
          solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_values);

          this->triangulation->execute_coarsening_and_refinement();

          this->dof_handler.distribute_dofs(*this->fe);
          this->locally_relevantRg_vector.reinit(this->dof_handler.n_dofs());

          solution_transfer.interpolate(old_Rg_values, this->locally_relevantRg_vector);

          this->construct_Rg_vector();*/

        } else {
          this->triangulation->execute_coarsening_and_refinement();
        }
      }

    }

    template <int dim>
    void  RefinementGlobal<dim>::print_convergence_table() const
    {
      this->convergence_table->omit_column_from_convergence_rate_evaluation("cycle");
      this->convergence_table->omit_column_from_convergence_rate_evaluation("cells");
      this->convergence_table->omit_column_from_convergence_rate_evaluation("DoFs");
      this->convergence_table->evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
      cout<<std::endl;
      this->convergence_table->write_text(std::cout);
      cout<<std::endl;

      std::ofstream tex_file(OUTPUT_PATH+"/convergence_table.tex");
      if (!tex_file.is_open()) {
        throw std::runtime_error("Failed to open convergence_table.tex for writing");
      }
      this->convergence_table->write_tex(tex_file);
      tex_file.close();
    }


    template <int dim>
    WeightedResidual<dim>::CellData::CellData(
      const FiniteElement<dim> &fe,
      const Quadrature<dim> &   quadrature,
      const Function<dim> &     right_hand_side,
      const MappingQ<dim> & mapping)
      : mapping(mapping),
        fe_values(mapping,
                  fe,
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
      : mapping(cell_data.mapping),
        fe_values(cell_data.mapping,
                  cell_data.fe_values.get_fe(),
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
      const Quadrature<dim - 1> &face_quadrature,
      const MappingQ<dim> & mapping)
      : mapping(mapping),
        fe_face_values_cell(mapping,
                            fe,
                            face_quadrature,
                            update_values | update_gradients |
                              update_JxW_values | update_normal_vectors)
      , fe_face_values_neighbor(mapping,
                                fe,
                                face_quadrature,
                                update_values | update_gradients |
                                  update_JxW_values | update_normal_vectors)
      , fe_subface_values_cell(mapping, fe, face_quadrature, update_gradients)
    {
      const unsigned int n_face_q_points = face_quadrature.size();

      jump_residual.resize(n_face_q_points);
      dual_weights.resize(n_face_q_points);
      cell_grads.resize(n_face_q_points);
      neighbor_grads.resize(n_face_q_points);
    }

    template <int dim>
    WeightedResidual<dim>::FaceData::FaceData(const FaceData &face_data)
      : mapping(face_data.mapping)
      , fe_face_values_cell(face_data.mapping,
                            face_data.fe_face_values_cell.get_fe(),
                            face_data.fe_face_values_cell.get_quadrature(),
                            update_values | update_gradients |
                              update_JxW_values | update_normal_vectors)
      , fe_face_values_neighbor(
          face_data.mapping,
          face_data.fe_face_values_neighbor.get_fe(),
          face_data.fe_face_values_neighbor.get_quadrature(),
          update_values | update_gradients | update_JxW_values |
            update_normal_vectors)
      , fe_subface_values_cell(
          face_data.mapping,
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
        const Vector<double> &     dual_weights,
        const MappingQ<dim> & mapping)
      : mapping(mapping),
        cell_data(primal_fe, primal_quadrature, rhs_function, mapping)
      , face_data(primal_fe, primal_face_quadrature, mapping)
      , primal_solution(primal_solution)
      , dual_weights(dual_weights)
    {}

    template <int dim>
    WeightedResidual<dim>::WeightedResidualScratchData::
      WeightedResidualScratchData(
        const WeightedResidualScratchData &scratch_data)
      : mapping(scratch_data.mapping)
      , cell_data(scratch_data.cell_data)
      , face_data(scratch_data.face_data)
      , primal_solution(scratch_data.primal_solution)
      , dual_weights(scratch_data.dual_weights)
    {}

    template <int dim>
    void output_cell_flags_to_vtk(const parallel::distributed::Triangulation<dim> &triangulation, unsigned int refinement_cycle, MPI_Comm & mpi_communicator) {
      Vector<float> cell_flags(triangulation.n_active_cells());
      unsigned int cell_index = 0;
      for (const auto &cell : triangulation.active_cell_iterators())
      {
        if (cell->is_locally_owned()) {
          if (cell->refine_flag_set())
            cell_flags[cell_index] = 1.0; // Refinement flag
          else if (cell->coarsen_flag_set())
            cell_flags[cell_index] = -1.0; // Coarsening flag
          else
            cell_flags[cell_index] = 0.0; // No flag
          ++cell_index;
        }
      }
      DataOut<dim> data_out;
      data_out.attach_triangulation(triangulation);
      data_out.add_data_vector(cell_flags, "CellFlags",DataOut<dim>::type_cell_data);
      data_out.build_patches();
      data_out.write_vtu_with_pvtu_record(
        OUTPUT_PATH+"/", "cell_flags", refinement_cycle, mpi_communicator, 2);

    }

    template <int dim>
    WeightedResidual<dim>::WeightedResidual(
      parallel::distributed::Triangulation<dim> &                           coarse_grid,
      const FiniteElement<dim> &                     primal_fe,
      const FiniteElement<dim> &                     dual_fe,
      const Quadrature<dim> &                        quadrature,
      const Quadrature<dim - 1> &                    face_quadrature,
      const Function<dim> &                          rhs_function,
      const Function<dim> &                          bv,
      const DualFunctional::DualFunctionalBase<dim> &dual_functional,
      const unsigned degree)
      : Base<dim>(coarse_grid)
      , PrimalSolver<dim>(coarse_grid,
                          primal_fe,
                          quadrature,
                          face_quadrature,
                          rhs_function,
                          bv,
                          degree)
      , DualSolver<dim>(coarse_grid,
                        dual_fe,
                        quadrature,
                        face_quadrature,
                        dual_functional,
                        rhs_function,
                        bv,
                        degree)
    {}

    template <int dim>
    void WeightedResidual<dim>::solve_problem()
    {
      solve_primal_problem();
      solve_dual_problem();
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
      TimerOutput::Scope t(this->computing_timer, "refine");

      // ERROR ESTIMATION
      Vector<float> error_indicators(this->triangulation->n_active_cells());
      error_indicators = 0.0;

      estimate_error(error_indicators);

      for (float &error_indicator : error_indicators)
        error_indicator = std::fabs(error_indicator);

      unsigned int nonzero_entries = 0;
      for (size_t i = 0; i<error_indicators.size(); ++i)
        if (std::fabs(error_indicators[i]) > 1e-8)
          nonzero_entries++;
      cout<<"Rank "<<this->this_mpi_process<<std::scientific<<" nonzero entries: "<<nonzero_entries<<std::endl;
      MPI_Barrier(this->mpi_communicator);

      // GRID REFINEMENT

      if (REFINEMENT_CRITERION == 3) {
        if (MANUAL_LIFTING_ON) {
          /*Vector<double> old_Rg_values = PrimalSolver<dim>::locally_relevant_Rg_vector;
          SolutionTransfer<dim> solution_transfer(PrimalSolver<dim>::dof_handler);
          solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_values);
          this->triangulation->refine_global(1);
          PrimalSolver<dim>::dof_handler.distribute_dofs(*(PrimalSolver<dim>::fe));
          PrimalSolver<dim>::locally_relevant_Rg_vector.reinit(PrimalSolver<dim>::dof_handler.n_dofs());
          solution_transfer.interpolate(old_Rg_values, PrimalSolver<dim>::locally_relevant_Rg_vector);
          PrimalSolver<dim>::construct_Rg_vector();
          DualSolver<dim>::locally_relevantRg_vector.reinit(DualSolver<dim>::dof_handler.n_dofs());*/
        } else {
          this->triangulation->refine_global(1);
        }
        return;
      }

      if (REFINEMENT_STRATEGY==1)
        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                          error_indicators,
                                                          TOP_FRACTION,
                                                          BOTTOM_FRACTION);
      else if (REFINEMENT_STRATEGY==2)
        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
                                                          error_indicators,
                                                          TOP_FRACTION,
                                                          BOTTOM_FRACTION);
      else if (REFINEMENT_STRATEGY==3) {
        // TODO: not avaliable in parallel
        /*GridRefinement::refine_and_coarsen_optimize(*this->triangulation,
                                                          error_indicators,
                                                          OPTIMIZE_ORDER);*/
      }
      else
        DEAL_II_NOT_IMPLEMENTED();

      if (MANUAL_LIFTING_ON) {    // TODO
        /*this->triangulation->prepare_coarsening_and_refinement();
        SolutionTransfer<dim> solution_transfer(PrimalSolver<dim>::dof_handler);

        Vector<double> old_Rg_values = PrimalSolver<dim>::locally_relevantRg_vector;
        solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_values);

        output_cell_flags_to_vtk(*this->triangulation, this->refinement_cycle, this->mpi_communicator);
        this->triangulation->execute_coarsening_and_refinement();

        PrimalSolver<dim>::dof_handler.distribute_dofs(*PrimalSolver<dim>::fe);
        PrimalSolver<dim>::locally_relevantRg_vector.reinit(PrimalSolver<dim>::dof_handler.n_dofs());

        solution_transfer.interpolate(old_Rg_values, PrimalSolver<dim>::locally_relevantRg_vector);

        PrimalSolver<dim>::construct_Rg_vector();

        DualSolver<dim>::locally_relevantRg_vector.reinit(DualSolver<dim>::dof_handler.n_dofs());*/
      } else {
        this->triangulation->prepare_coarsening_and_refinement();
        this->triangulation->execute_coarsening_and_refinement();
      }


    }

    template <int dim>
    void WeightedResidual<dim>::output_solution()
    {
      TimerOutput::Scope t(this->computing_timer, "output");

      PETScWrappers::MPI::Vector completely_distributed_dual_solution_on_primal(PrimalSolver<dim>::locally_owned_dofs,
                                                                                 PrimalSolver<dim>::mpi_communicator); // Constructor for NO GHOST elements
      FETools::interpolate(DualSolver<dim>::dof_handler,
                           DualSolver<dim>::locally_relevant_solution,
                           PrimalSolver<dim>::dof_handler,
                           PrimalSolver<dim>::linear_system_ptr->hanging_node_constraints,
                           completely_distributed_dual_solution_on_primal);

      DataOut<dim> data_out;
      data_out.attach_dof_handler(PrimalSolver<dim>::dof_handler);
      data_out.add_data_vector(PrimalSolver<dim>::locally_relevant_solution, "uh",DataOut<dim, dim>::type_dof_data);

      if (MANUAL_LIFTING_ON) { //TODO
        //data_out.add_data_vector(PrimalSolver<dim>::homogeneous_solution, "uh0",DataOut<dim, dim>::type_dof_data);
        //data_out.add_data_vector(PrimalSolver<dim>::locally_relevantRg_vector, "Rg",DataOut<dim, dim>::type_dof_data);
      }

      data_out.add_data_vector(completely_distributed_dual_solution_on_primal, "zh",DataOut<dim, dim>::type_dof_data);

      Vector<double> boundary_ids(this->triangulation->n_active_cells());
      for (const auto &cell : this->triangulation->active_cell_iterators())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            if (cell->is_locally_owned())
              boundary_ids[cell->active_cell_index()] = cell->face(face)->boundary_id();
      data_out.add_data_vector(boundary_ids, "boundary_ids",DataOut<dim, dim>::type_cell_data);

      Vector<float> subdomain(this->triangulation->n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = this->triangulation->locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");


      data_out.build_patches(PrimalSolver<dim>::mapping, 1, DataOut<dim,dim>::CurvedCellRegion::curved_inner_cells);

      const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record(
        OUTPUT_PATH+"/", "solution", this->refinement_cycle, this->mpi_communicator, 2);


      /*GridOut grid_out;
      GridOutFlags::Msh msh_flags(true, true);
      grid_out.set_flags(msh_flags);
      grid_out.write_msh(*this->triangulation, OUTPUT_PATH+"/final_mesh.msh");*/

    }

    template <int dim>
    void WeightedResidual<dim>::estimate_error(Vector<float> &error_indicators) const {

      PETScWrappers::MPI::Vector completely_distributed_primal_solution_on_dual(
        DualSolver<dim>::locally_owned_dofs,
        DualSolver<dim>::mpi_communicator);   // non-ghosted

      FETools::interpolate(PrimalSolver<dim>::dof_handler,
                           PrimalSolver<dim>::locally_relevant_solution, // ghosted
                           DualSolver<dim>::dof_handler,
                           DualSolver<dim>::linear_system_ptr->hanging_node_constraints,
                           completely_distributed_primal_solution_on_dual);  // non-ghosted

      PETScWrappers::MPI::Vector locally_relevant_primal_solution_on_dual(
        DualSolver<dim>::locally_owned_dofs,
        DualSolver<dim>::locally_relevant_dofs,
        DualSolver<dim>::mpi_communicator);   // ghosted
      locally_relevant_primal_solution_on_dual = completely_distributed_primal_solution_on_dual;

      // --------

      PETScWrappers::MPI::Vector completely_distributed_dual_weights(
        DualSolver<dim>::locally_owned_dofs,
        DualSolver<dim>::mpi_communicator);  // non-ghosted

      FETools::interpolation_difference(
        DualSolver<dim>::dof_handler,
        DualSolver<dim>::linear_system_ptr->hanging_node_constraints,
        DualSolver<dim>::locally_relevant_solution,    // ghosted
        PrimalSolver<dim>::dof_handler,
        PrimalSolver<dim>::linear_system_ptr->hanging_node_constraints,
        completely_distributed_dual_weights);   // non-ghosted

      PETScWrappers::MPI::Vector locally_relevant_dual_weights(
        DualSolver<dim>::locally_owned_dofs,
        DualSolver<dim>::locally_relevant_dofs,
        DualSolver<dim>::mpi_communicator);   // ghosted
      locally_relevant_dual_weights = completely_distributed_dual_weights;

      // --------

      auto & primal_solution = locally_relevant_primal_solution_on_dual;
      auto & dual_weights = locally_relevant_dual_weights;


      // ------------------------------------------------------------
      // LOCAL ESTIMATE: integrate over cells
      // ------------------------------------------------------------

      FEValues<dim> fe_values(DualSolver<dim>::mapping, *DualSolver<dim>::fe,
                              *DualSolver<dim>::quadrature,
                              update_values | update_gradients | update_quadrature_points | update_JxW_values);

      const unsigned int n_q_points = DualSolver<dim>::quadrature->size();

      std::vector<Tensor<1,dim>> cell_Rg_plus_uh0hat_gradients(n_q_points);
      std::vector<Tensor<1,dim>> cell_dual_weights_gradients(n_q_points);
      std::vector<double> cell_rhs_values(n_q_points);
      std::vector<double> cell_dual_weights(n_q_points);

      double sum;

      for (const auto &cell : DualSolver<dim>::dof_handler.active_cell_iterators()){
        if (cell->is_locally_owned()) {
          fe_values.reinit(cell);
          fe_values.get_function_gradients(primal_solution, cell_Rg_plus_uh0hat_gradients);
          fe_values.get_function_gradients(dual_weights, cell_dual_weights_gradients);
          fe_values.get_function_values(dual_weights, cell_dual_weights);
          PrimalSolver<dim>::rhs_function->value_list(fe_values.get_quadrature_points(), cell_rhs_values);

          // Numerically approximate the integral of the scalar product between the gradients of the two
          sum = 0;

          for (unsigned int p = 0; p < n_q_points; ++p) {
            sum -=  eps_r * eps_0 *
                    ((cell_Rg_plus_uh0hat_gradients[p] * cell_dual_weights_gradients[p]  )   // Scalar product btw Tensors
                      * fe_values.JxW(p));
            sum += (cell_rhs_values[p] * cell_dual_weights[p] * fe_values.JxW(p));
          }
          error_indicators(cell->active_cell_index()) += sum;
        }
      }

      // BEGIN Output
      DataOut<dim> data_out;
      data_out.attach_dof_handler(DualSolver<dim>::dof_handler);
      data_out.add_data_vector(error_indicators, "error_indicators", DataOut<dim, dim>::type_cell_data);
      data_out.add_data_vector(primal_solution, "primal_solution", DataOut<dim, dim>::type_dof_data);
      data_out.add_data_vector(dual_weights, "dual_weights", DataOut<dim, dim>::type_dof_data);
      data_out.build_patches();
      const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record(
        "./../results/", "error_indicators", this->refinement_cycle, this->mpi_communicator, 4);
      // END Output

      this->pcout<<"error_indicators.size() = "<<error_indicators.size()<<std::endl;
      unsigned int nonzero_entries = 0;
      for (size_t i = 0; i<error_indicators.size(); ++i)
        if (std::fabs(error_indicators[i]) > 1e-8)
          nonzero_entries++;
      cout<<"Rank "<<this->this_mpi_process<<std::scientific<<" nonzero entries: "<<nonzero_entries<<std::endl;
      MPI_Barrier(this->mpi_communicator);

      //double local_error = std::accumulate(error_indicators.begin(), error_indicators.end(), 0.0);
      // Accumulate only locally owned cells' contributions
      double local_error = 0.0;
      for (const auto &cell : DualSolver<dim>::dof_handler.active_cell_iterators()) {
        if (cell->is_locally_owned()) {
          local_error += error_indicators(cell->active_cell_index());
        }
      }

      cout<<"Rank "<<this->this_mpi_process<<std::scientific<<" local error: "<<std::fabs(local_error)<<std::endl;

      double global_error = 0.0;
      MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, this->mpi_communicator);
      cout<<std::scientific<<"Common error (MPI_Allreduce): "<<std::fabs(global_error)<<std::endl;
      MPI_Barrier(this->mpi_communicator);

      double estimated_error = Utilities::MPI::sum(local_error, this->mpi_communicator);

      PrimalSolver<dim>::convergence_table->add_value("est err",std::fabs(estimated_error));
      PrimalSolver<dim>::convergence_table->set_scientific("est err",true);

      CSVLogger::getInstance().addColumn("est err", to_string_with_precision(std::fabs(estimated_error),15));

    }

    /*template <int dim>
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
                                    dual_weights,
                                    DualSolver<dim>::mapping),
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

      double estimated_error = std::accumulate(error_indicators.begin(),
                                   error_indicators.end(),
                                   0.);


      PrimalSolver<dim>::convergence_table->add_value("est err",std::fabs(estimated_error));
      PrimalSolver<dim>::convergence_table->set_scientific("est err",true);

      CSVLogger::getInstance().addColumn("est err", to_string_with_precision(std::fabs(estimated_error),15));
    }*/

    template <int dim>
    void WeightedResidual<dim>::update_convergence_table() {
      if (this->this_mpi_process == 0)
      {
        PrimalSolver<dim>::convergence_table->add_value("cycle", this->refinement_cycle);
        PrimalSolver<dim>::convergence_table->add_value("cells", this->triangulation->n_active_cells());
        PrimalSolver<dim>::convergence_table->add_value("DoFs", PrimalSolver<dim>::dof_handler.n_dofs());

        CSVLogger& logger = CSVLogger::getInstance();
        logger.addColumn("cycle", std::to_string(this->refinement_cycle));
        logger.addColumn("cells", std::to_string(this->triangulation->n_active_cells()));
        logger.addColumn("DoFs", std::to_string(PrimalSolver<dim>::dof_handler.n_dofs()));
      }
    }

    template <int dim>
    void WeightedResidual<dim>::print_convergence_table() const
    {
      // No convergence rates. Makes no sense
      cout<<std::endl;
      PrimalSolver<dim>::convergence_table->write_text(std::cout);
      cout<<std::endl;
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
        sum += ((cell_data.rhs_values[p] + eps_0 * eps_r * cell_data.cell_laplacians[p]) *
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
        face_integral += eps_0 * eps_r *
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
            face_integral += eps_0 * eps_r *
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
    template class WeightedResidual<2>;
  } // namespace LaplaceSolver
} // namespace IonPropulsion