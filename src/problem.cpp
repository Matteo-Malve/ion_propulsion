#include "problem.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

using namespace dealii;
using std::cout;
using std::endl;

template <int dim>
Problem<dim>::Problem() : primal_dof_handler(triangulation), 
                          dual_dof_handler(triangulation),
                          primal_fe(1), 
                          dual_fe(1), 
                          cycle(0) {}

template <int dim>
void Problem<dim>::run() {
	create_mesh();
	while (cycle <= NUM_REFINEMENT_CYCLES) {
		cout << endl << "Cycle " << cycle << ':' << endl;
		std::cout << "   Number of active cells:       "<< triangulation.n_active_cells() << std::endl;
		// Primal --------------
		cout<<"   Primal:"<<endl;
		setup_primal_system();
		assemble_primal_system();
		solve_primal();
		output_primal_results();

		// --------------------  
    if(cycle==NUM_REFINEMENT_CYCLES){
      // Last cycle primal -------------------- 
			cout << "   Primal " << " [FINAL]" << ':' << endl;
      setup_primal_system();
      assemble_primal_system();
      solve_primal();
      output_primal_results();
    } else {
      // Dual ---------------    
      cout<<"   Dual:"<<endl;
      setup_dual_system();
      assemble_dual_system();
      solve_dual();
      output_dual_results();

			// Error estimation and grid refinement ---------------  
      cout<<"   Error estimation and Mesh refinement:"<<endl;
      refine_mesh();
    }
    ++cycle;
	}
}

template <int dim>
void Problem<dim>::create_mesh() {
	cout << endl << "Reading file: " << extract_mesh_name() << endl;
	std::ifstream input_file(PATH_TO_MESH);
	GridIn<2>       grid_in;
	grid_in.attach_triangulation(triangulation);
	grid_in.read_msh(input_file);

  for (unsigned int i = 0; i < NUM_PRELIMINARY_REF; ++i) {
		Vector<float> criteria(triangulation.n_active_cells());
		cout  << "Active cells " << triangulation.n_active_cells() << endl;
		unsigned int ctr = 0;

    // Threshold
    const double max_thickness = 2. * L;
    const double min_thickness = 1.05 * L;
    const double D = min_thickness + (max_thickness-min_thickness)/(NUM_PRELIMINARY_REF-1)*(NUM_PRELIMINARY_REF-1-i);
    cout<<"D = "<<D<<endl;

		for (auto &cell : triangulation.active_cell_iterators()) {                
			const Point<dim> c = cell->center();
        if(c[1]<D && std::abs(c[0])<D)
          criteria[ctr++] = 1;
        else
          criteria[ctr++] = 0;
		}
		GridRefinement::refine(triangulation, criteria, 0.5);
		triangulation.execute_coarsening_and_refinement();
		cout<<"   Executed preliminary coarsening and refinement"<<endl;
	}
}

// -----------------------------------------
// PRIMAL
// -----------------------------------------

template <int dim>
void Problem<dim>::setup_primal_system() {
  primal_dof_handler.distribute_dofs(primal_fe);
  
	//Rg_vector.reinit(primal_dof_handler.n_dofs());
	uh0.reinit(primal_dof_handler.n_dofs());
  primal_solution.reinit(primal_dof_handler.n_dofs());
  primal_rhs.reinit(primal_dof_handler.n_dofs());

  if (cycle == 0){
    // Initialize Rg_dof_values for the first time
    Rg_dof_values.reinit(primal_dof_handler.n_dofs());

    // Interpolate boundary values only once
    std::map<types::global_dof_index, double> emitter_boundary_values;
    VectorTools::interpolate_boundary_values(primal_dof_handler, 1, Functions::ConstantFunction<dim>(20000.), emitter_boundary_values);

    for (const auto &boundary_value : emitter_boundary_values)
      Rg_dof_values(boundary_value.first) = boundary_value.second;
  }

  primal_constraints.clear();
	DoFTools::make_hanging_node_constraints(primal_dof_handler, primal_constraints);
	primal_constraints.close();

  DynamicSparsityPattern dsp(primal_dof_handler.n_dofs(),primal_dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(primal_dof_handler,
                                  dsp
                                  //,primal_constraints,
                                  ///*keep_constrained_dofs = */ false
                                  );
  primal_constraints.condense(dsp);
  primal_sparsity_pattern.copy_from(dsp);
  primal_system_matrix.reinit(primal_sparsity_pattern);

  std::cout << "      Number of degrees of freedom: " << primal_dof_handler.n_dofs()<<std::endl;
  
}

template <int dim>
void Problem<dim>::assemble_primal_system() {
  const QGauss <dim> quadrature(primal_fe.degree + 1);
  FEValues<dim> fe_values(primal_fe,
                          quadrature,
                          update_values | update_gradients | update_quadrature_points |
                          update_JxW_values);

  const unsigned int dofs_per_cell = primal_fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const unsigned int n_q_points = quadrature.size();
  std::vector<Tensor<1, dim>> rg_gradients(n_q_points);

  for (const auto &cell : primal_dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;

    // Compute gradients of Rg at quadrature points
    fe_values.get_function_gradients(Rg_dof_values, rg_gradients);

    // Compute A_loc and rhs_loc
    for (const unsigned int q_index : fe_values.quadrature_point_indices()){ 
      for (const unsigned int i : fe_values.dof_indices()){
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) += eps_r*eps_0*
                              (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                                fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                                fe_values.JxW(q_index));           // dx

        cell_rhs(i) -= eps_r * eps_0 *
                        (fe_values.shape_grad(i, q_index) *   // grad phi_i(x_q)
                          rg_gradients[q_index] *             // grad_Rg(x_q)
                          fe_values.JxW(q_index));            // dx
        
      }
    }

    // Local to global
    cell->get_dof_indices(local_dof_indices);

    // Chatgpt suggests to hide loop:
    primal_rhs.add(local_dof_indices, cell_rhs);

    //primal_constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, primal_system_matrix, primal_rhs);
    for (const unsigned int i : fe_values.dof_indices()){
      for (const unsigned int j : fe_values.dof_indices())
        primal_system_matrix.add(local_dof_indices[i], 
                          local_dof_indices[j],
                          cell_matrix(i, j)); 
    }           
  }

  // Apply boundary values
  std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;
  // Emitter
  VectorTools::interpolate_boundary_values(primal_dof_handler,1, Functions::ZeroFunction<dim>(), emitter_boundary_values);
  // Collector
  VectorTools::interpolate_boundary_values(primal_dof_handler,2, Functions::ZeroFunction<dim>(), collector_boundary_values);
  
  // Condense constraints
  primal_constraints.condense(primal_system_matrix);
  primal_constraints.condense(primal_rhs);

  MatrixTools::apply_boundary_values(emitter_boundary_values, primal_system_matrix, primal_solution, primal_rhs);
  MatrixTools::apply_boundary_values(collector_boundary_values, primal_system_matrix, primal_solution, primal_rhs);

}

template <int dim>
void Problem<dim>::solve_primal() {
  // Solver setup
  const unsigned int it_max = 1e+4;
  const double rel_tol = 1.e-6*primal_rhs.l2_norm();
  const double abs_tol = 1.e-12;
  const double tol = abs_tol + rel_tol;
  SolverControl            solver_control(it_max, tol);
  SolverCG<Vector<double>> solver(solver_control);

  // Solve linear system --> uh0
  double relaxation_parameter = 1.2;
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(primal_system_matrix, relaxation_parameter);
  solver.solve(primal_system_matrix, uh0, primal_rhs, preconditioner);

  // uh [primal_solution]
  primal_solution = uh0;
  primal_solution += Rg_dof_values;

  // Distribute constraints
	primal_constraints.distribute(uh0);
  primal_constraints.distribute(primal_solution);

  cout<<"      Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

template <int dim>
void Problem<dim>::output_primal_results() {
	ElectricFieldPostprocessor<dim> electric_field_postprocessor;
  HomogeneousFieldPostprocessor<dim> homogeneous_field_postprocessor;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(primal_dof_handler);
  data_out.add_data_vector(primal_solution, "Potential");
	data_out.add_data_vector(uh0, "uh0");
	data_out.add_data_vector(Rg_dof_values, "Rg");
  data_out.add_data_vector(primal_solution, electric_field_postprocessor);
  data_out.add_data_vector(uh0, homogeneous_field_postprocessor);
  
  data_out.build_patches(); // mapping

  std::string filename;
  std::string meshName = extract_mesh_name();
	filename = std::string("Rg_manual") + "-" + "primal" + "-" + meshName + "-" + Utilities::int_to_string(cycle, 1) + ".vtk";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtk(output);
}

// -----------------------------------------
// DUAL
// -----------------------------------------

template <int dim>
void Problem<dim>::setup_dual_system() {
  dual_dof_handler.distribute_dofs(dual_fe);

  dual_solution.reinit(dual_dof_handler.n_dofs());
  dual_rhs.reinit(dual_dof_handler.n_dofs());

  if (cycle == 0){
    // Initialize Rg_dof_values for the first time
    Rg_dual_dof_values.reinit(dual_dof_handler.n_dofs());

    // Interpolate boundary values only once
    std::map<types::global_dof_index, double> emitter_boundary_values;
    VectorTools::interpolate_boundary_values(dual_dof_handler, 1, Functions::ConstantFunction<dim>(20000.), emitter_boundary_values);

    for (const auto &boundary_value : emitter_boundary_values)
      Rg_dual_dof_values(boundary_value.first) = boundary_value.second;
  }

  dual_constraints.clear();
	DoFTools::make_hanging_node_constraints(dual_dof_handler, dual_constraints);                   
	dual_constraints.close();

  DynamicSparsityPattern dsp(dual_dof_handler.n_dofs(),dual_dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dual_dof_handler, dsp);
  dual_constraints.condense(dsp);

  dual_sparsity_pattern.copy_from(dsp);
  dual_system_matrix.reinit(dual_sparsity_pattern);

  std::cout << "      Number of degrees of freedom: " << dual_dof_handler.n_dofs()<< std::endl;
}

template <int dim>
void Problem<dim>::assemble_dual_system() {
  // ---------------------------------------
  // RHS
  // ---------------------------------------
  Point<dim> evaluation_point;
  if (NUM_PRELIMINARY_REF == 4)
    evaluation_point = Point<dim>(0.00025, 0.0005);  
  else{
    cout<<"Choose an evaluation point suitable for "<<NUM_PRELIMINARY_REF<<"initial refinements"<<endl;
    abort();
  }
  //PointValueEvaluation<dim> dual_functional(evaluation_point);
  BoundaryFluxEvaluation<dim> dual_functional(1);  // Pass boundary ID, e.g., 1
  //FaceBoundaryFluxEvaluation<dim> dual_functional(1);  // Pass boundary ID, e.g., 1
  dual_functional.assemble_rhs(dual_dof_handler, dual_rhs);

  // ---------------------------------------
  // MATRIX
  // ---------------------------------------
  const QGauss<dim> quadrature(dual_dof_handler.get_fe().degree + 1);
  FEValues<dim> fe_values(dual_fe,
                          quadrature,
                          update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = dual_fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dual_dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_matrix = 0;

    // Compute A_loc
    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) += eps_r*eps_0*
                  (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                    fe_values.JxW(q_index));           // dx
    
    // Local to global
    cell->get_dof_indices(local_dof_indices);
    for (const unsigned int i : fe_values.dof_indices())
      for (const unsigned int j : fe_values.dof_indices())
        dual_system_matrix.add(local_dof_indices[i], 
                          local_dof_indices[j],
                          cell_matrix(i, j)); 
  }

  // Apply boundary values
  std::map<types::global_dof_index, double> emitter_and_collector_boundary_values;
  VectorTools::interpolate_boundary_values(dual_dof_handler,1, Functions::ZeroFunction<dim>(), emitter_and_collector_boundary_values);
  VectorTools::interpolate_boundary_values(dual_dof_handler,2, Functions::ZeroFunction<dim>(), emitter_and_collector_boundary_values);
  
  // Condense constraints
  dual_constraints.condense(dual_system_matrix);
  dual_constraints.condense(dual_rhs);

  MatrixTools::apply_boundary_values(emitter_and_collector_boundary_values, dual_system_matrix, dual_solution, dual_rhs);
}

template <int dim>
void Problem<dim>::solve_dual() {
  const unsigned int it_max = 1e+4;
  const double rel_tol = 1.e-6*dual_rhs.l2_norm();
  const double abs_tol = 1.e-12;

  const double tol = abs_tol + rel_tol;
  SolverControl            solver_control(it_max, tol);
  SolverCG<Vector<double>> solver(solver_control);

  double relaxation_parameter = 1.2;
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(dual_system_matrix, relaxation_parameter);

  solver.solve(dual_system_matrix, dual_solution, dual_rhs, preconditioner);

  dual_constraints.distribute(dual_solution);

  cout<<"      Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

template <int dim>
void Problem<dim>::output_dual_results() {
	ElectricFieldPostprocessor<dim> electric_field_postprocessor;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dual_dof_handler);
  data_out.add_data_vector(dual_solution, "Potential");
  data_out.add_data_vector(dual_solution, electric_field_postprocessor);
  Vector<double> constraint_indicator(dual_dof_handler.n_dofs());
  for (unsigned int i = 0; i < dual_dof_handler.n_dofs(); ++i)
    constraint_indicator(i) = dual_constraints.is_constrained(i) ? 1.0 : 0.0;

  data_out.add_data_vector(constraint_indicator, "constraints");
  data_out.build_patches(5); // mapping

  std::string filename;
  std::string meshName = extract_mesh_name();

	filename =  std::string("Rg_manual") + "-" + "dual-" + meshName + "-" + Utilities::int_to_string(cycle, 1) + ".vtk";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtk(output);
}

// -----------------------------------------
// GOAL ORIENTED REFINEMENT
// -----------------------------------------

template <int dim>
double Problem<dim>::estimate_error(Vector<float> &error_indicators) const{
	// ------------------------------------------------------------      
  // PROJECTIONS: for both LOCAL and GLOBAL estimates

  Vector<double> primal_homogeneous_solution_on_dual_space(dual_dof_handler.n_dofs());
  const auto sol_size = primal_homogeneous_solution_on_dual_space.size();
  FETools::interpolate(primal_dof_handler,
                        uh0, 
                        dual_dof_handler, 
                        dual_constraints,
                        primal_homogeneous_solution_on_dual_space);

  dual_constraints.distribute(primal_homogeneous_solution_on_dual_space);

  // Subtract from dual solution its projection on the primal solution FE space
  Vector<double> dual_weights(dual_dof_handler.n_dofs());
  FETools::interpolation_difference(dual_dof_handler,
                                    dual_constraints,
                                    dual_solution,
                                    primal_dof_handler,
                                    primal_constraints,
                                    dual_weights);               // zh-∏zh
  dual_constraints.distribute(dual_weights); 

  // ------------------------------------------------------------      
  // RETRIEVE LIFTING: Rg + uh0hat
  // ! Here we were summing up Rg and then forgetting the vector forever

  // uh0^hat + Rg
  Vector<double> Rg_plus_uh0hat(dual_dof_handler.n_dofs());
  Rg_plus_uh0hat = primal_homogeneous_solution_on_dual_space;
  Rg_plus_uh0hat += Rg_dual_dof_values;
  dual_constraints.distribute(Rg_plus_uh0hat);

  // ------------------------------------------------------------      
  // LOCAL ESTIMATE: integrate over cells
  const QGauss<dim> quadrature(dual_dof_handler.get_fe().degree + 1);
  FEValues<dim> fe_values(dual_fe,
                          quadrature,
                          update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = quadrature.size();

  std::vector<Tensor<1,dim>> cell_primal_gradients(n_q_points);
  std::vector<Tensor<1,dim>> cell_dual_gradients(n_q_points);

  double sum;

  for (const auto &cell : dual_dof_handler.active_cell_iterators()){

    fe_values.reinit(cell);
    fe_values.get_function_gradients(Rg_plus_uh0hat, cell_primal_gradients);
    fe_values.get_function_gradients(dual_weights, cell_dual_gradients);

    // Numerically approximate the integral of the scalar product between the gradients of the two
    sum = 0;
    
    for (unsigned int p = 0; p < n_q_points; ++p) {
        sum +=  eps_r * eps_0 *                        
                ((cell_primal_gradients[p] * cell_dual_gradients[p]  )   // Scalar product btw Tensors
                  * fe_values.JxW(p));
    }
    error_indicators(cell->active_cell_index()) -= sum;

  }

  // ------------------------------------------------------------      
  // GLOBAL ESTIMATE
  
  double dx = 0.0; double lx = 0.0;

  // Compute   u' A z              NOTE: z is actually the dual weights (zh-∏zh)
  Vector<double> temp(dual_dof_handler.n_dofs());
  dual_system_matrix.vmult(temp,dual_weights);    // A z      NB: A is condensed
  dual_constraints.distribute(temp);

  for(unsigned int i=0;i<sol_size;++i)
      dx+=primal_homogeneous_solution_on_dual_space(i)*temp(i);           // u' A z

  // Compute   a(Rg,φ)
  const unsigned int dofs_per_cell = dual_fe.n_dofs_per_cell();

  Vector<double> F(dual_dof_handler.n_dofs());
  Vector<double> cell_F(dofs_per_cell);
  std::vector <types::global_dof_index> local_dof_indices(dofs_per_cell);

  // fe_values declared above for local estimate
  std::vector<Tensor<1, dim>> rg_gradients(n_q_points);

  for (const auto &cell: dual_dof_handler.active_cell_iterators()) {
    cell_F = 0;
    fe_values.reinit(cell);

    // Compute gradients of Rg at quadrature points
    fe_values.get_function_gradients(Rg_dual_dof_values, rg_gradients);

    // assemble a(Rg,φ)
    for (const unsigned int q_index : fe_values.quadrature_point_indices()){ 
      for (const unsigned int i : fe_values.dof_indices())
        cell_F(i) += eps_r * eps_0*
                        (fe_values.shape_grad(i, q_index) *     // grad phi_i(x_q)
                        rg_gradients[q_index] *                               // grad_Rg(x_q)
                        fe_values.JxW(q_index));                // dx
    }
    // Local to Global
    cell->get_dof_indices(local_dof_indices);
    F.add(local_dof_indices, cell_F);
  }
  dual_constraints.condense(F);
  dual_constraints.distribute(F);

  // Compute     z' F    or    z•F
  for(unsigned int i=0;i<dual_weights.size();i++)
    lx += dual_weights(i) * F(i);

  // Return             η = |r(z)|  = | - z' F - u' A z |
  // or, precisely,     η = |r(zk-∏zk)|  = | - (zk-∏zk)' F - uh0' A (zk-∏zk) |
  return std::abs(-lx -dx);
}

template <int dim>
void Problem<dim>::refine_mesh() {
  // Prepare vector to store error values
  Vector<float> error_indicators(triangulation.n_active_cells());
  // Compute local errors and global error estimate
  double global_error = estimate_error(error_indicators);
  // Sum contribution of each cell's local error to get a global estimate
  double global_error_as_sum_of_cell_errors=0.0;
  for(unsigned int i=0; i<error_indicators.size(); i++)
      global_error_as_sum_of_cell_errors += error_indicators[i];
  global_error_as_sum_of_cell_errors = std::abs(global_error_as_sum_of_cell_errors);

  // Output the two derived global estimates
  cout<<"      Global error = " <<  global_error << endl
      <<"      Global error as sum of cells' errors = " << global_error_as_sum_of_cell_errors << endl;

  // Take absolute value of each error
  for (float &error_indicator : error_indicators) {
      error_indicator = std::fabs(error_indicator);
  }
  
	GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,error_indicators, 0.8, 0); // CHANGED FROM: 0.8, 0.02
  
  // Prepare the solution transfer object
  SolutionTransfer<dim> solution_transfer(primal_dof_handler);
  
  SolutionTransfer<dim> dual_solution_transfer(dual_dof_handler);
  // take a copy of the solution vector
  Vector<double> old_Rg_dof_values(Rg_dof_values);
  Vector<double> old_dual_Rg_dof_values(Rg_dual_dof_values); 
  // Prepare for refinement (older versions of deal.II)
  solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_dof_values);
  dual_solution_transfer.prepare_for_coarsening_and_refinement(old_dual_Rg_dof_values);
  // Perform the refinement
  triangulation.execute_coarsening_and_refinement();
  // Reinitialize the DoFHandler for the refined mesh
  primal_dof_handler.distribute_dofs(primal_fe);
  dual_dof_handler.distribute_dofs(dual_fe);

  
  // Reinitialize Rg_dof_values to match the new DoF layout after refinement
  Rg_dof_values.reinit(primal_dof_handler.n_dofs());
  Rg_dual_dof_values.reinit(dual_dof_handler.n_dofs());
  // Transfer the old values to the new DoFs, accounting for hanging nodes
  solution_transfer.interpolate(old_Rg_dof_values, Rg_dof_values);
  dual_solution_transfer.interpolate(old_dual_Rg_dof_values, Rg_dual_dof_values);

  // Handle boundary conditions again (for hanging nodes)
  std::map<types::global_dof_index, double> emitter_boundary_values;
  VectorTools::interpolate_boundary_values(primal_dof_handler, 1, Functions::ConstantFunction<dim>(20000.), emitter_boundary_values);
  for (const auto &boundary_value : emitter_boundary_values)
    Rg_dof_values(boundary_value.first) = boundary_value.second;
  std::map<types::global_dof_index, double> dual_emitter_boundary_values;
  VectorTools::interpolate_boundary_values(dual_dof_handler, 1, Functions::ConstantFunction<dim>(20000.), dual_emitter_boundary_values);
  for (const auto &boundary_value : dual_emitter_boundary_values)
    Rg_dual_dof_values(boundary_value.first) = boundary_value.second;
}




// #######################################
// Template initialization
// #######################################
template class Problem<2>;