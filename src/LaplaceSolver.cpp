#include "LaplaceSolver.h"

#include <Constants.h>

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

    template <int dim>
    void Base<dim>::checkpoint()
    {
      std::cout << "--- Writing checkpoint... ---" << std::endl << std::endl;

      /*{
        std::ofstream checkpoint_file(OUTPUT_PATH+"/"+"tmp.checkpoint_ion_propulsion");
        AssertThrow(checkpoint_file,
                    ExcMessage(
                      "Could not write to the <tmp.checkpoint_ion_propulsion> file."));

        boost::archive::text_oarchive archive(checkpoint_file);

        archive << *this;
      }*/

      //particle_handler.prepare_for_serialization();
      //triangulation->save(OUTPUT_PATH+"/"+"tmp.checkpoint");
      std::ofstream                 checkpoint_file(OUTPUT_PATH+"/"+"checkpoint-mesh-"+std::to_string(this->refinement_cycle));
      boost::archive::text_oarchive archive(checkpoint_file);
      triangulation->save(archive,0);
      /*
      std::list<std::string> tmp_checkpoint_files;
      for (const auto &dir_entry : std::filesystem::directory_iterator(OUTPUT_PATH))
        if (dir_entry.is_regular_file() &&
            (dir_entry.path().filename().string().find("tmp.checkpoint") == 0))
          tmp_checkpoint_files.push_back(dir_entry.path().filename().string());

      for (const std::string &filename : tmp_checkpoint_files) {
        std::string new_name = filename.substr(4, std::string::npos);
        std::filesystem::rename(OUTPUT_PATH + "/" + filename, OUTPUT_PATH + "/" + new_name);
      }*/
    }

    template <int dim>
    void Solver<dim>::restart()
      {
        /*{
          std::ifstream checkpoint_file(OUTPUT_PATH+"/"+"checkpoint_ion_propulsion");
          AssertThrow(checkpoint_file,
                      ExcMessage(
                        "Could not read from the <checkpoint_ion_propulsion> file."));

          boost::archive::text_iarchive archive(checkpoint_file);
          archive >> *this;
        }

      this->triangulation->load(OUTPUT_PATH+"/"+"checkpoint");
      //particle_handler.deserialize();

      dof_handler.reinit(*this->triangulation);

      GridOut grid_out;
      GridOutFlags::Msh msh_flags(true, true);
      grid_out.set_flags(msh_flags);
      grid_out.write_msh(*this->triangulation, OUTPUT_PATH+"/re-imported_mesh.msh");*/

      cout<<"Enter restart() function"<<std::endl;
      std::ifstream checkpoint_file("../checkpoint-mesh");
      AssertThrow(checkpoint_file,
                  ExcMessage(
                    "Could not read from the <checkpoint-mesh> file."));
      cout<<"Checkpoint file found"<<std::endl;
      boost::archive::text_iarchive archive(checkpoint_file);

      cout<<"Archive built. Ready to load traingulation:"<<std::endl;
      this->triangulation->load(archive,1);
      cout<<"Mesh loaded."<<std::endl;

      dof_handler.reinit(*this->triangulation);
    }


    // ------------------------------------------------------
    // Solver
    // ------------------------------------------------------

    template <int dim>
    Solver<dim>::Solver(Triangulation<dim> &       triangulation,
                        const FiniteElement<dim> & fe,
                        const Quadrature<dim> &    quadrature,
                        const Quadrature<dim - 1> &face_quadrature,
                        const Function<dim> &      boundary_values,
                        const unsigned degree)
      : Base<dim>(triangulation)
      , fe(&fe)
      , quadrature(&quadrature)
      , face_quadrature(&face_quadrature)
      , dof_handler(triangulation)
      , boundary_values(&boundary_values)
      , mapping(degree)
    {}


    template <int dim>
    Solver<dim>::~Solver()
    {
      dof_handler.clear();
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
    void Solver<dim>::solve_problem()
    {
      dof_handler.distribute_dofs(*fe);
      homogeneous_solution.reinit(dof_handler.n_dofs());
      solution.reinit(dof_handler.n_dofs());

      if(MANUAL_LIFTING_ON) {
        if (this->refinement_cycle == 0) {
          Rg_vector.reinit(dof_handler.n_dofs());
          construct_Rg_vector();
        }
      }


      // Solve homogeneous system with lifting in the rhs
      linear_system_ptr = std::make_unique<LinearSystem>(dof_handler);
      assemble_linear_system(*linear_system_ptr);

      linear_system_ptr->solve(homogeneous_solution);
      solution = homogeneous_solution;

      if (LOAD_FROM_SETUP != 0 && LOAD_FROM_SETUP != 11)
        compute_second_order_flux(*linear_system_ptr);

      // Retrieve lifting
      if(MANUAL_LIFTING_ON)
        retrieve_Rg();


    }


    template <int dim>
    void Solver<dim>::postprocess(
      const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      std::pair<std::string, double> possible_pair = postprocessor(dof_handler, solution,*this->triangulation);
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
                      AssemblyScratchData(*fe, *quadrature,mapping),
                      AssemblyCopyData());
      //linear_system.hanging_node_constraints.condense(linear_system.matrix);

      std::map<types::global_dof_index, double> boundary_value_map;
      interpolate_boundary_values(boundary_value_map);

      rhs_task.join();
      linear_system.hanging_node_constraints.condense(linear_system.rhs);

      MatrixTools::apply_boundary_values(boundary_value_map,
                                         linear_system.matrix,
                                         homogeneous_solution,
                                         linear_system.rhs);

    }

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
              (scratch_data.fe_values.shape_grad(i, q_point) *
               scratch_data.fe_values.shape_grad(j, q_point) *
               scratch_data.fe_values.JxW(q_point));

      cell->get_dof_indices(copy_data.local_dof_indices);
    }



    template <int dim>
    void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data,
                                           LinearSystem &linear_system) const
    {
      linear_system.hanging_node_constraints.distribute_local_to_global(
        copy_data.cell_matrix,
        copy_data.local_dof_indices,
        linear_system.matrix);

      for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i)
        for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j) {
          /*linear_system.matrix.add(copy_data.local_dof_indices[i],
                                   copy_data.local_dof_indices[j],
                                   copy_data.cell_matrix(i, j));*/
          linear_system.Umatrix.add(copy_data.local_dof_indices[i],
                                   copy_data.local_dof_indices[j],
                                   copy_data.cell_matrix(i, j));
        }

    }

    template <int dim>
    Solver<dim>::AssemblyScratchData::AssemblyScratchData(
      const FiniteElement<dim> &fe,
      const Quadrature<dim> &   quadrature,
      MappingQ<dim> & mapping)
      : mapping(mapping),
        fe_values(mapping, fe, quadrature, update_gradients | update_JxW_values)
    {}


    template <int dim>
    Solver<dim>::AssemblyScratchData::AssemblyScratchData(
      const AssemblyScratchData &scratch_data)
      : mapping(scratch_data.mapping),
        fe_values(scratch_data.mapping,
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  update_gradients | update_JxW_values)
    {}


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
      Umatrix.reinit(sparsity_pattern);
      rhs.reinit(dof_handler.n_dofs());
    }



    template <int dim>
    void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const
    {
      SolverControl            solver_control(5000, 1e-12);
      SolverCG<Vector<double>> cg(solver_control);

      PreconditionSSOR<SparseMatrix<double>> preconditioner;
      preconditioner.initialize(matrix, 1.2);

      cg.solve(matrix, solution, rhs, preconditioner);
      cout<<"Checkpoint"<<std::endl;
      hanging_node_constraints.distribute(solution);

      cout<<"Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<std::endl;
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
                                    const Function<dim> &      boundary_values,
                                    const unsigned degree)
      : Base<dim>(triangulation)
      , Solver<dim>(triangulation,
                    fe,
                    quadrature,
                    face_quadrature,
                    boundary_values,
                    degree)
      , rhs_function(&rhs_function)
    {}



    template <int dim>
    void PrimalSolver<dim>::output_solution() const
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler(this->dof_handler);
      data_out.add_data_vector(this->solution, "uh",DataOut<dim, dim>::type_dof_data);

      if(MANUAL_LIFTING_ON) {
        data_out.add_data_vector(this->homogeneous_solution, "uh0",DataOut<dim, dim>::type_dof_data);
        data_out.add_data_vector(this->Rg_vector, "Rg",DataOut<dim, dim>::type_dof_data);
      }

      Vector<double> rhs_function_values(this->dof_handler.n_dofs());
      VectorTools::interpolate(this->dof_handler, *this->rhs_function, rhs_function_values);
      data_out.add_data_vector(rhs_function_values, "rhs_function",DataOut<dim, dim>::type_dof_data);

      Vector<double> uex_function_values(this->dof_handler.n_dofs());
      VectorTools::interpolate(this->dof_handler, *this->boundary_values, uex_function_values);
      data_out.add_data_vector(uex_function_values, "u_ex",DataOut<dim, dim>::type_dof_data);

      Vector<double> boundary_ids(this->triangulation->n_active_cells());
      for (const auto &cell : this->triangulation->active_cell_iterators())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            boundary_ids[cell->active_cell_index()] = cell->face(face)->boundary_id();
      data_out.add_data_vector(boundary_ids, "boundary_ids",DataOut<dim, dim>::type_cell_data);

      data_out.build_patches(this->mapping, 1,DataOut<dim,dim>::CurvedCellRegion::curved_inner_cells);

      std::ofstream out(OUTPUT_PATH+"/"+"solution-" + std::to_string(this->refinement_cycle) +
                        ".vtu");

      data_out.write(out, DataOutBase::vtu);

    }



    template <int dim>
    void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const
    {
      FEValues<dim> fe_values(this->mapping,
                              *this->fe,
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

        rhs_function->value_list(fe_values.get_quadrature_points(),
                                 rhs_values);
        if(MANUAL_LIFTING_ON)
          fe_values.get_function_gradients(this->Rg_vector, rg_gradients);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                            rhs_values[q_point] *               // f((x_q)
                            fe_values.JxW(q_point));            // dx
            if(MANUAL_LIFTING_ON) {
              cell_rhs(i) -= eps_r * eps_0 *
                      (fe_values.shape_grad(i, q_point) *   // grad phi_i(x_q)
                        rg_gradients[q_point] *             // grad_Rg(x_q)
                        fe_values.JxW(q_point));            // dx
            }
          }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          rhs(local_dof_indices[i]) += cell_rhs(i);
      }
    }

    template <int dim>
    void PrimalSolver<dim>::compute_second_order_flux(typename Solver<dim>::LinearSystem &linear_system) {
      // Extract set of DoFs on the emitter
      IndexSet complete_index_set = this->dof_handler.locally_owned_dofs();
      auto e_index_set = DoFTools::extract_boundary_dofs(this->dof_handler,
                                      ComponentMask(),
                                      std::set<types::boundary_id>({1}));
      double flux = 0.;
      //Vector<double> Au(this->dof_handler.n_dofs());
      Au.reinit(this->dof_handler.n_dofs());

      // Compute flux = - sum((Au-b)(e_nodes))
      // Step 1: Au-b
      if (MANUAL_LIFTING_ON) {
        linear_system.Umatrix.vmult(Au, this->solution);
        Au-=linear_system.rhs;
        Vector<double> ARg(this->dof_handler.n_dofs());
        linear_system.Umatrix.vmult(ARg, this->Rg_vector);
        Au+=ARg;
      } else { // TODO: for now works only on zero Dirichlet BCs
        linear_system.Umatrix.vmult(Au, this->solution);
        Au-=linear_system.rhs;
      }

      // Step 2: for i in e_nodes DO flux -= (Au-b)[i]
      for (auto index = e_index_set.begin(); index != e_index_set.end(); ++index) {
        flux += Au(*index);
      }
      std::cout << std::scientific << std::setprecision(12)
                  << "   Cons. flux = " << flux << std::endl;

      this->conservative_flux = flux;


    }

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
      if (LOAD_FROM_SETUP==0) {
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 1, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 2, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 3, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
      } else if (LOAD_FROM_SETUP==11) {
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 1, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 2, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 3, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 4, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
      } else {
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 1, *(this->boundary_values), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 2, *(this->boundary_values), boundary_value_map);
        VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler, 9, *(this->boundary_values), boundary_value_map); //TODO: remove 9
      }
      for (const auto &boundary_value : boundary_value_map)
        this->Rg_vector(boundary_value.first) = boundary_value.second;

      //cout<< std::scientific << std::setprecision(12)<<"Exact point value at EVALUATION POINT : "<<this->boundary_values->value(Point<dim>(EVALUATION_POINT_X,EVALUATION_POINT_Y))<<std::endl;

      side_task.join();
      hanging_node_constraints.close();
      hanging_node_constraints.distribute(this->Rg_vector);
    }

    template <int dim>
    void PrimalSolver<dim>::interpolate_boundary_values(std::map<types::global_dof_index, double> & boundary_value_map) {
      if(MANUAL_LIFTING_ON) {
        if (LOAD_FROM_SETUP==0) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
        } else if (LOAD_FROM_SETUP == 11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ZeroFunction<dim>(),boundary_value_map);
        } else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,Functions::ZeroFunction<dim>(),boundary_value_map);
        }
      } else {
        if (LOAD_FROM_SETUP==0) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ConstantFunction<dim>(Ve),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ConstantFunction<dim>(Vc),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ConstantFunction<dim>(Vc),boundary_value_map);
        } else if (LOAD_FROM_SETUP==11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ConstantFunction<dim>(Ve),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ConstantFunction<dim>(Ve),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ConstantFunction<dim>(Vc),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ConstantFunction<dim>(Vc),boundary_value_map);
        }
        else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,*this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,*this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,*this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,*this->boundary_values,boundary_value_map);
        }
      }

    }

    // ------------------------------------------------------
    // DualSolver
    // ------------------------------------------------------

    template <int dim>
    void assemble_conservative_flux_rhs(
      const DoFHandler<dim> &dof_handler,
      Vector<double> &rhs,
      SparseMatrix<double> & Umatrix,
      const FiniteElement<dim> & fe,
      const Quadrature<dim> & quadrature,
      const Function<dim> & special_rhs_function,
      const Function<dim> & special_boundary_values,
      AffineConstraints<double> & hanging_node_constraints

      ) {
      // ------------------------------------------------------
      // Primal-like Rg
      // ------------------------------------------------------
      Vector<double> special_Rg_vector(dof_handler.n_dofs());

      std::map<types::global_dof_index, double> boundary_value_map;
      if (LOAD_FROM_SETUP==0) {
        VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 3, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
      } else if (LOAD_FROM_SETUP==11) {
        VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ConstantFunction<dim>(Ve), boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 3, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 4, Functions::ConstantFunction<dim>(Vc), boundary_value_map);
      } else {
        VectorTools::interpolate_boundary_values(dof_handler, 1, special_boundary_values, boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 2, special_boundary_values, boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 3, special_boundary_values, boundary_value_map);
        VectorTools::interpolate_boundary_values(dof_handler, 9, special_boundary_values, boundary_value_map);
      }
      for (const auto &boundary_value : boundary_value_map)
        special_Rg_vector(boundary_value.first) = boundary_value.second;
      hanging_node_constraints.distribute(special_Rg_vector);

      // ------------------------------------------------------
      // Primal-like b
      // ------------------------------------------------------
      Vector<double> b(dof_handler.n_dofs());
      FEValues<dim> fe_values(fe,
                              quadrature,
                              update_values | update_gradients | update_quadrature_points |
                              update_JxW_values);
      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
      const unsigned int n_q_points    = quadrature.size();
      Vector<double>                       cell_rhs(dofs_per_cell);
      std::vector<double>                  rhs_values(n_q_points);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Tensor<1, dim>>          rg_gradients(n_q_points);

      for (const auto &cell : dof_handler.active_cell_iterators()){
        cell_rhs = 0;
        fe_values.reinit(cell);
        special_rhs_function.value_list(fe_values.get_quadrature_points(),
                                 rhs_values);
        if(MANUAL_LIFTING_ON)
          fe_values.get_function_gradients(special_Rg_vector, rg_gradients);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
          for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                            rhs_values[q_point] *               // f((x_q)
                            fe_values.JxW(q_point));            // dx
            if(MANUAL_LIFTING_ON) {
              cell_rhs(i) -= eps_r * eps_0 *
                      (fe_values.shape_grad(i, q_point) *   // grad phi_i(x_q)
                        rg_gradients[q_point] *             // grad_Rg(x_q)
                        fe_values.JxW(q_point));            // dx
            }
          }
        }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          b(local_dof_indices[i]) += cell_rhs(i);
      }

      // Check //TODO: DELETE
      unsigned int b_nonzero_elements = 0;
      for (size_t i = 0; i<rhs.size(); ++i)
        if (std::fabs(b(i))>1.e-10)
          b_nonzero_elements++;
      cout<<"b's nonzero elements: "<<b_nonzero_elements<<std::endl;


      // ------------------------------------------------------
      // Put together for conservative
      // ------------------------------------------------------
      IndexSet complete_index_set = dof_handler.locally_owned_dofs();
      auto e_index_set = DoFTools::extract_boundary_dofs(dof_handler,
                                      ComponentMask(),
                                      std::set<types::boundary_id>({1}));

      Vector<double> v(dof_handler.n_dofs());
      Vector<double> A_col_i(dof_handler.n_dofs());
      unsigned int nonzero_elements = 0;

      for (size_t i=0; i<b.size(); ++i) {
        v.reinit(dof_handler.n_dofs());
        A_col_i.reinit(dof_handler.n_dofs());
        v(i)=1.;
        Umatrix.vmult(v, A_col_i);

        Vector<double> ARg(dof_handler.n_dofs());// JUST ADDED
        Umatrix.vmult(ARg, special_Rg_vector);
        A_col_i+=ARg;

        //A_col_i-=b;

        //rhs(i)=A_col_i(i)-b(i);
        for (auto index = e_index_set.begin(); index != e_index_set.end(); ++index) {
          rhs(i) += A_col_i(*index);
        }
        //rhs(i) = std::accumulate(A_col_i.begin(),A_col_i.end(), 0.);

        if (std::fabs(rhs(i))>1.e-10)
          nonzero_elements++;
      }
      cout<<"Nonzero elements (rhs perspective): "<<nonzero_elements<<std::endl;

    }

    template <int dim>
    DualSolver<dim>::DualSolver(
      Triangulation<dim> &                           triangulation,
      const FiniteElement<dim> &                     fe,
      const Quadrature<dim> &                        quadrature,
      const Quadrature<dim - 1> &                    face_quadrature,
      const DualFunctional::DualFunctionalBase<dim> &dual_functional,
      const Function<dim> &                          special_rhs_function,
      const Function<dim> &                          special_boundary_values,
      const unsigned degree)
      : Base<dim>(triangulation)
      , Solver<dim>(triangulation,
                    fe,
                    quadrature,
                    face_quadrature,
                    boundary_values,
                    degree)
      , special_rhs_function(&special_rhs_function)
      , special_boundary_values(&special_boundary_values)
      , dual_functional(&dual_functional)
    {}

    template <int dim>
    void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const
    {
      if (DUAL_FUNCTIONAL==3)
        assemble_conservative_flux_rhs(
          this->dof_handler,
          rhs,
          this->linear_system_ptr->Umatrix,
          *this->fe,
          *this->quadrature,
          *special_rhs_function,
          *special_boundary_values,
           this->linear_system_ptr->hanging_node_constraints);
      else
        dual_functional->assemble_rhs(this->dof_handler, rhs);
      //this->conservative_flux_rhs(rhs);
    }

    template <int dim>
    void DualSolver<dim>::interpolate_boundary_values(std::map<types::global_dof_index, double> & boundary_value_map) {
      if(MANUAL_LIFTING_ON) {
        if (LOAD_FROM_SETUP==0) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
        } else if (LOAD_FROM_SETUP == 11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ZeroFunction<dim>(),boundary_value_map);
        } else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,Functions::ZeroFunction<dim>(),boundary_value_map);
        }
      } else {
        if (LOAD_FROM_SETUP==0) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
        } else if (LOAD_FROM_SETUP==11) {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,Functions::ZeroFunction<dim>(),boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,4,Functions::ZeroFunction<dim>(),boundary_value_map);
        }
        else {
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,1,this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,2,this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,3,this->boundary_values,boundary_value_map);
          VectorTools::interpolate_boundary_values(this->mapping,this->dof_handler,9,this->boundary_values,boundary_value_map);
        }
      }

    }






    // Template instantiation
    template class Base<2>;
    template class Solver<2>;
    template class PrimalSolver<2>;
    template class DualSolver<2>;

  } // namespace LaplaceSolver
} // namespace IonPropulsion