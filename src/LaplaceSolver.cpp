#include "LaplaceSolver.h"
#include "Constants.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{

    // ------------------------------------------------------
    // Base
    // ------------------------------------------------------
    template <int dim>
    Base<dim>::Base(Triangulation<dim> &coarse_grid)
      : triangulation(&coarse_grid)
      , convergence_table(std::make_shared<ConvergenceTable>())
      , refinement_cycle(numbers::invalid_unsigned_int)
    {}

    template <int dim>
    void Base<dim>::set_refinement_cycle(const unsigned int cycle)
    {
      refinement_cycle = cycle;
    }

    // ------------------------------------------------------
    // Solver
    // ------------------------------------------------------

    template <int dim>
    Solver<dim>::Solver(Triangulation<dim> &       triangulation,
                        const FiniteElement<dim> & fe,
                        const Quadrature<dim> &    quadrature,
                        const Quadrature<dim - 1> &face_quadrature,
                        const Function<dim> &      boundary_values)
      : Base<dim>(triangulation)
      , fe(&fe)
      , quadrature(&quadrature)
      , face_quadrature(&face_quadrature)
      , dof_handler(triangulation)
      , boundary_values(&boundary_values)
    {}


    template <int dim>
    Solver<dim>::~Solver()
    {
      dof_handler.clear();
    }


    template <int dim>
    void Solver<dim>::solve_problem()
    {
      dof_handler.distribute_dofs(*fe);
      homogeneous_solution.reinit(dof_handler.n_dofs());
      solution.reinit(dof_handler.n_dofs());

      if (this->refinement_cycle == 0) {
        Rg_vector.reinit(dof_handler.n_dofs());
        construct_Rg_vector();
      }

      linear_system_ptr = std::make_unique<LinearSystem>(dof_handler);
      assemble_linear_system(*linear_system_ptr);
      linear_system_ptr->solve(homogeneous_solution /*solution*/);

      // Retrieve lifting
      solution = homogeneous_solution;
      retrieve_Rg();

      // Compute flux
      //compute_flux();

    }

    template <int dim>
    void Solver<dim>::compute_flux() {

      std::unique_ptr<const Function<dim>> rhs_function = std::make_unique<Functions::ConstantFunction<dim>>(1.);
      // TODO: Generalize this.
      // --> rhs_function available only in PrimalSolver

      //SparseMatrix<double>      flux_matrix;
      //Vector<double>            flux_rhs;
      double flux = 0.;

      const QGauss<dim-1> face_quadrature(dof_handler.get_fe().degree + 1);
      FEFaceValues<dim> fe_face_values(*this->fe,
                              face_quadrature,
                              update_values | update_gradients | update_quadrature_points |
                              update_JxW_values);
      const unsigned int dofs_per_face = this->fe->n_dofs_per_face(); // TODO: devo reinizializzare? Cambai con hanging nodes?
      const unsigned int n_face_q_points = face_quadrature.size();

      std::vector<double>                  rhs_values(n_face_q_points);
      std::vector<Tensor<1, dim>>          rg_gradients(n_face_q_points);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_face);

      // Cycle on all cells
      for (const auto &cell : dof_handler.active_cell_iterators()) {

        // Determine if cell is on the EMITTER
        for (const auto &face : cell->face_iterators())
          if (face->at_boundary() && face->boundary_id() == 1) {

            fe_face_values.reinit(cell,face);

            // Evaluate: gradRg & F
            fe_face_values.get_function_gradients(this->Rg_vector, rg_gradients);
            rhs_function->value_list(fe_face_values.get_quadrature_points(),rhs_values);


            for (const unsigned int q_index : fe_face_values.quadrature_point_indices()) {
              for (const unsigned int i : fe_face_values.dof_indices()){
                for (const unsigned int j : fe_face_values.dof_indices())
                  flux +=
                    eps_r * eps_0 *
                    solution(local_dof_indices[j]) *             // w_j
                    (fe_face_values.shape_grad(i, q_index) *   // grad phi_i(x_q)
                    fe_face_values.shape_grad(j, q_index) *  // grad phi_j(x_q)
                    fe_face_values.JxW(q_index));                     // dx
                // NOTE: Inverted signs for rhs contributions
                flux -=
                  (fe_face_values.shape_value(i, q_index) *   // phi_i(x_q)
                  rhs_values[q_index] *                         // f(x_q)
                  fe_face_values.JxW(q_index));                      //dx
                /*flux +=
                  eps_r * eps_0 *
                  (fe_face_values.shape_grad(i, q_index) *   // grad phi_i(x_q)
                  rg_gradients[q_index] *                      // grad_Rg(x_q)
                  fe_face_values.JxW(q_index));                     // dx*/
            }
          }

        } // end if boundary
      } // end cell loop


      // Output the result
      std::cout << std::scientific << std::setprecision(12)
                << "   2.Flux=" << flux << std::endl;


      /*

      // Apply boundary values
      std::map<types::global_dof_index, double> boundary_values;
      VectorTools::interpolate_boundary_values(primal_dof_handler,1, Functions::ZeroFunction<dim>(), boundary_values);
      VectorTools::interpolate_boundary_values(primal_dof_handler,9, Functions::ZeroFunction<dim>(), boundary_values);
      //VectorTools::interpolate_boundary_values(primal_dof_handler,1, exact_solution_function, boundary_values);
      //VectorTools::interpolate_boundary_values(primal_dof_handler,9, exact_solution_function, boundary_values);

      // Condense constraints
      primal_constraints.condense(primal_system_matrix);
      primal_constraints.condense(primal_rhs);

      MatrixTools::apply_boundary_values(boundary_values, primal_system_matrix, uh0, primal_rhs); */
    }

    template <int dim>
    void Solver<dim>::update_convergence_table() {
      this->convergence_table->add_value("cycle", this->refinement_cycle);
      this->convergence_table->add_value("cells", this->triangulation->n_active_cells());
      this->convergence_table->add_value("DoFs", this->dof_handler.n_dofs());

      CSVLogger& logger = CSVLogger::getInstance();
      logger.addColumn("cycle", std::to_string(this->refinement_cycle));
      logger.addColumn("cells", std::to_string(this->triangulation->n_active_cells()));
      logger.addColumn("DoFs", std::to_string(this->dof_handler.n_dofs()));
    }


    template <int dim>
    void Solver<dim>::postprocess(
      const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      Evaluation::ExtraData extra_data;
      if (linear_system_ptr) {
        extra_data.matrix = &linear_system_ptr->matrix;
        extra_data.rhs = &linear_system_ptr->rhs;
      }

      std::pair<std::string, double> possible_pair = postprocessor(dof_handler, solution, &extra_data);

      if (possible_pair.first != "null") {
        Base<dim>::convergence_table->add_value(possible_pair.first, possible_pair.second);
        Base<dim>::convergence_table->set_scientific(possible_pair.first, true);
        CSVLogger::getInstance().addColumn(possible_pair.first, to_string_with_precision(possible_pair.second,15));
      }


    }


    template <int dim>
    unsigned int Solver<dim>::n_dofs() const
    {
      return dof_handler.n_dofs();
    }


    // The following few functions and constructors are verbatim
    // copies taken from step-13:
    template <int dim>
    void Solver<dim>::assemble_linear_system(LinearSystem &linear_system)
    {
      Threads::Task<void> rhs_task =
        Threads::new_task(&Solver<dim>::assemble_rhs, *this, linear_system.rhs);

      auto worker =
        [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
               AssemblyScratchData &scratch_data,
               AssemblyCopyData &   copy_data) {
          this->local_assemble_matrix(cell, scratch_data, copy_data);
        };

      auto copier = [this, &linear_system](const AssemblyCopyData &copy_data) {
        this->copy_local_to_global(copy_data, linear_system);
      };

      WorkStream::run(dof_handler.begin_active(),
                      dof_handler.end(),
                      worker,
                      copier,
                      AssemblyScratchData(*fe, *quadrature),
                      AssemblyCopyData());

      linear_system.hanging_node_constraints.condense(linear_system.matrix);

      std::map<types::global_dof_index, double> boundary_value_map;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),  /*Now assigned by Rg in construction*/
                                               boundary_value_map);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(),  /*Now assigned by Rg in construction*/
                                               boundary_value_map);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               9,
                                               Functions::ZeroFunction<dim>(),  /*Now assigned by Rg in construction*/
                                               boundary_value_map);
      rhs_task.join();
      linear_system.hanging_node_constraints.condense(linear_system.rhs);
      MatrixTools::apply_boundary_values(boundary_value_map,
                                         linear_system.matrix,
                                         homogeneous_solution /*solution*/,
                                         linear_system.rhs);
    }


    template <int dim>
    Solver<dim>::AssemblyScratchData::AssemblyScratchData(
      const FiniteElement<dim> &fe,
      const Quadrature<dim> &   quadrature)
      : fe_values(fe, quadrature, update_gradients | update_JxW_values)
    {}


    template <int dim>
    Solver<dim>::AssemblyScratchData::AssemblyScratchData(
      const AssemblyScratchData &scratch_data)
      : fe_values(scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  update_gradients | update_JxW_values)
    {}


    template <int dim>
    void Solver<dim>::local_assemble_matrix(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      AssemblyScratchData &                                 scratch_data,
      AssemblyCopyData &                                    copy_data) const
    {
      const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
      const unsigned int n_q_points    = quadrature->size();

      copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);

      copy_data.local_dof_indices.resize(dofs_per_cell);

      scratch_data.fe_values.reinit(cell);

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            copy_data.cell_matrix(i, j) += eps_r * eps_0 *
              (scratch_data.fe_values.shape_grad(i, q_point) *          // grad phi_i(x_q)
               scratch_data.fe_values.shape_grad(j, q_point) *        // grad phi_j(x_q)
               scratch_data.fe_values.JxW(q_point));                    // dx

      cell->get_dof_indices(copy_data.local_dof_indices);
    }



    template <int dim>
    void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data,
                                           LinearSystem &linear_system) const
    {
      for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i)
        for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j)
          linear_system.matrix.add(copy_data.local_dof_indices[i],
                                   copy_data.local_dof_indices[j],
                                   copy_data.cell_matrix(i, j));
    }


      template <int dim>
      Solver<dim>::LinearSystem::LinearSystem(const DoFHandler<dim> &dof_handler)
      {
        hanging_node_constraints.clear();

        void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) =
          &DoFTools::make_hanging_node_constraints;

        // Start a side task then continue on the main thread
        Threads::Task<void> side_task =
          Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints);

        DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
        DoFTools::make_sparsity_pattern(dof_handler, dsp);



        // Wait for the side task to be done before going further
        side_task.join();

        hanging_node_constraints.close();
        hanging_node_constraints.condense(dsp);
        sparsity_pattern.copy_from(dsp);

        matrix.reinit(sparsity_pattern);
        rhs.reinit(dof_handler.n_dofs());
      }



      template <int dim>
      void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const   // Note: We will pass homogeneous_solution
      {
        SolverControl            solver_control(1000, 1e-12);
        SolverCG<Vector<double>> cg(solver_control);

        PreconditionSSOR<SparseMatrix<double>> preconditioner;
        preconditioner.initialize(matrix, 1.2);

        cg.solve(matrix, solution, rhs, preconditioner);

        hanging_node_constraints.distribute(solution);

      //cout<<"Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<std::endl;
      }

    // ------------------------------------------------------
    // PrimalSolver
    // ------------------------------------------------------

    template <int dim>
    PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &       triangulation,
                                    const FiniteElement<dim> & fe,
                                    const Quadrature<dim> &    quadrature,
                                    const Quadrature<dim - 1> &face_quadrature,
                                    const Function<dim> &      rhs_function,
                                    const Function<dim> &      boundary_values)
      : Base<dim>(triangulation)
      , Solver<dim>(triangulation,
                    fe,
                    quadrature,
                    face_quadrature,
                    boundary_values)
      , rhs_function(&rhs_function)
    {}

    template <int dim>
    void PrimalSolver<dim>::construct_Rg_vector() {
      AffineConstraints<double> hanging_node_constraints;
      hanging_node_constraints.clear();
      void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) =
          &DoFTools::make_hanging_node_constraints;
      // Start a side task then continue on the main thread
      Threads::Task<void> side_task =
        Threads::new_task(mhnc_p, this->dof_handler, hanging_node_constraints);

      std::map<types::global_dof_index, double> boundary_value_map;
      VectorTools::interpolate_boundary_values(this->dof_handler, 0, *(this->boundary_values), boundary_value_map);
      VectorTools::interpolate_boundary_values(this->dof_handler, 1, *(this->boundary_values), boundary_value_map);
      VectorTools::interpolate_boundary_values(this->dof_handler, 9, *(this->boundary_values), boundary_value_map);
      for (const auto &boundary_value : boundary_value_map)
        this->Rg_vector(boundary_value.first) = boundary_value.second;

      side_task.join();
      hanging_node_constraints.close();
      hanging_node_constraints.distribute(this->Rg_vector);
    }



    template <int dim>
    void PrimalSolver<dim>::output_solution() const
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler(this->dof_handler);
      data_out.add_data_vector(this->homogeneous_solution, "uh0");
      data_out.add_data_vector(this->solution, "uh");
      data_out.add_data_vector(this->Rg_vector, "Rg");

      Vector<double> boundary_ids(this->triangulation->n_active_cells());
      for (const auto &cell : this->triangulation->active_cell_iterators())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            boundary_ids[cell->active_cell_index()] = cell->face(face)->boundary_id();
      data_out.add_data_vector(boundary_ids, "boundary_ids");

      data_out.build_patches();

      std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                        ".vtu");
      data_out.write(out, DataOutBase::vtu);

    }



    template <int dim>
    void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const
    {
      FEValues<dim> fe_values(*this->fe,
                              *this->quadrature,
                              update_values | update_gradients | update_quadrature_points |
                                update_JxW_values);
      const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
      const unsigned int n_q_points    = this->quadrature->size();

      Vector<double>                       cell_rhs(dofs_per_cell);
      std::vector<double>                  rhs_values(n_q_points);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Tensor<1, dim>>          rg_gradients(n_q_points);

      for (const auto &cell : this->dof_handler.active_cell_iterators())
        {
          cell_rhs = 0;

          fe_values.reinit(cell);
          fe_values.get_function_gradients(this->Rg_vector, rg_gradients);
          rhs_function->value_list(fe_values.get_quadrature_points(),
                                   rhs_values);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                              rhs_values[q_point] *               // f((x_q)
                              fe_values.JxW(q_point));            // dx
              cell_rhs(i) -=   eps_r * eps_0 *
                        (fe_values.shape_grad(i, q_point) *   // grad phi_i(x_q)
                          rg_gradients[q_point] *             // grad_Rg(x_q)
                          fe_values.JxW(q_point));            // dx
            }
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    // ------------------------------------------------------
    // DualSolver
    // ------------------------------------------------------
    template <int dim>
    DualSolver<dim>::DualSolver(
      Triangulation<dim> &                           triangulation,
      const FiniteElement<dim> &                     fe,
      const Quadrature<dim> &                        quadrature,
      const Quadrature<dim - 1> &                    face_quadrature,
      const DualFunctional::DualFunctionalBase<dim> &dual_functional)
      : Base<dim>(triangulation)
      , Solver<dim>(triangulation,
                    fe,
                    quadrature,
                    face_quadrature,
                    boundary_values)
      , dual_functional(&dual_functional)
    {}



    template <int dim>
    void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const
    {
      dual_functional->assemble_rhs(this->dof_handler, rhs);
    }

    // Template instantiation
    template class Base<2>;
    template class Solver<2>;
    template class PrimalSolver<2>;
    template class DualSolver<2>;

  } // namespace LaplaceSolver
} // namespace IonPropulsion