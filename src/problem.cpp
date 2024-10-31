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
                          dual_fe(2)
                          {}

template <int dim>
void Problem<dim>::run() {
  // Import mesh, perform preliminary refinements
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
    test_convergence();
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
      // Error evaluation and grid refinement --------------- 
      cout<<"   Error estimation:"<<endl;
      estimate_error();
      output_dual_results();
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
  /*
  for (unsigned int i = 0; i < NUM_PRELIMINARY_REF; ++i) {
		Vector<float> criteria(triangulation.n_active_cells());
		cout  << "Active cells " << triangulation.n_active_cells() << endl;
		unsigned int ctr = 0;

    // Threshold
    const double max_thickness = 2. * l;
    const double min_thickness = 1.05 * l;
    const double D = min_thickness + (max_thickness-min_thickness)/(NUM_PRELIMINARY_REF-1)*(NUM_PRELIMINARY_REF-1-i);
    cout<<"D = "<<D<<endl;

		for (auto &cell : triangulation.active_cell_iterators()) {                
			const Point<dim> c = cell->center();
        if(std::abs(c[1])<D && std::abs(c[0])<D)
          criteria[ctr++] = 1;
        else
          criteria[ctr++] = 0;
		}
		GridRefinement::refine(triangulation, criteria, 0.5);
		triangulation.execute_coarsening_and_refinement();
		cout<<"   Executed preliminary coarsening and refinement"<<endl;
	}*/
  //triangulation.refine_global(2);
}

// -----------------------------------------
// PRIMAL
// -----------------------------------------

template <int dim>
void Problem<dim>::setup_primal_system() {
  primal_dof_handler.distribute_dofs(primal_fe);
  
	uh0.reinit(primal_dof_handler.n_dofs());
  uh.reinit(primal_dof_handler.n_dofs());
  primal_rhs.reinit(primal_dof_handler.n_dofs());

  if (cycle == 0){
    // Initialize Rg_primal for the first time
    Rg_primal.reinit(primal_dof_handler.n_dofs());

    // Interpolate boundary values only once
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(primal_dof_handler, 1, exact_solution_function, boundary_values);
    VectorTools::interpolate_boundary_values(primal_dof_handler, 9, exact_solution_function, boundary_values);

    for (const auto &boundary_value : boundary_values)
      Rg_primal(boundary_value.first) = boundary_value.second;
   
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
  const QGauss <dim> quadrature(4);
  FEValues<dim> fe_values(primal_fe,
                          quadrature,
                          update_values | update_gradients | update_quadrature_points |
                          update_JxW_values);

  const unsigned int dofs_per_cell = primal_fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  //std::vector<double> rhs_values(quadrature.size());
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const unsigned int n_q_points = quadrature.size();
  std::vector<Tensor<1, dim>> rg_gradients(n_q_points);

  for (const auto &cell : primal_dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;

    //rhs_function.value_list(fe_values.get_quadrature_points(),rhs_values);

    // Compute gradients of Rg at quadrature points
    fe_values.get_function_gradients(Rg_primal, rg_gradients);

    // Compute A_loc and rhs_loc
    for (const unsigned int q_index : fe_values.quadrature_point_indices()){ 
      for (const unsigned int i : fe_values.dof_indices()){
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) += eps_r * eps_0 *
                              (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                                fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                                fe_values.JxW(q_index));           // dx
        const auto &x_q = fe_values.quadrature_point(q_index);
        cell_rhs(i) += (fe_values.shape_value(i, q_index) *   // phi_i(x_q)
                          rhs_function.value(x_q) *           // f(x_q)
                          fe_values.JxW(q_index));            //dx
        cell_rhs(i) -= eps_r * eps_0 *
                        (fe_values.shape_grad(i, q_index) *   // grad phi_i(x_q)
                          rg_gradients[q_index] *             // grad_Rg(x_q)
                          fe_values.JxW(q_index));            // dx
      }
    }

    // Local to global
    cell->get_dof_indices(local_dof_indices);

    //primal_constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, primal_system_matrix, primal_rhs);
    for (const unsigned int i : fe_values.dof_indices()){
      primal_rhs(local_dof_indices[i]) += cell_rhs(i);
      for (const unsigned int j : fe_values.dof_indices())
        primal_system_matrix.add(local_dof_indices[i], 
                          local_dof_indices[j],
                          cell_matrix(i, j)); 
    }           
  }

  // Apply boundary values
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(primal_dof_handler,1, Functions::ZeroFunction<dim>(), boundary_values);
  VectorTools::interpolate_boundary_values(primal_dof_handler,9, Functions::ZeroFunction<dim>(), boundary_values);
  // floor gets Homogeneous_Neumann

  // Condense constraints
  primal_constraints.condense(primal_system_matrix);
  primal_constraints.condense(primal_rhs);

  MatrixTools::apply_boundary_values(boundary_values, primal_system_matrix, uh0, primal_rhs); 
  cout<<"      Applied BCs"<<endl;
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

  // uh [uh]
  uh = uh0;
  uh += Rg_primal;

  // Distribute constraints
	primal_constraints.distribute(uh0);
  primal_constraints.distribute(uh);

  cout<<"      Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

template <int dim>
void Problem<dim>::output_primal_results() {
	ElectricFieldPostprocessor<dim> electric_field_postprocessor;
  HomogeneousFieldPostprocessor<dim> homogeneous_field_postprocessor;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(primal_dof_handler);
  data_out.add_data_vector(uh, "Potential");
	data_out.add_data_vector(uh0, "uh0");
	data_out.add_data_vector(Rg_primal, "Rg");
  data_out.add_data_vector(uh, electric_field_postprocessor);
  data_out.add_data_vector(uh0, homogeneous_field_postprocessor);

  Vector<double> boundary_ids(triangulation.n_active_cells());
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->at_boundary())
        boundary_ids[cell->active_cell_index()] = cell->face(face)->boundary_id();
  data_out.add_data_vector(boundary_ids, "boundary_ids");
  
  Vector<double> u_ex(primal_dof_handler.n_dofs());
  VectorTools::project(primal_dof_handler, 
                      primal_constraints, 
                      QGauss<dim>(7),  // Quadrature rule (degree + 1 for accuracy)
                      exact_solution_function,          // Analytical function Rg
                      u_ex);
  data_out.add_data_vector(u_ex, "u_ex");

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

  zh.reinit(dual_dof_handler.n_dofs());
  dual_rhs.reinit(dual_dof_handler.n_dofs());

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
  PointValueEvaluation<dim> dual_functional(evaluation_point);
  //BoundaryFluxEvaluation<dim> dual_functional(1);  // Pass boundary ID, e.g., 1
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

  MatrixTools::apply_boundary_values(emitter_and_collector_boundary_values, dual_system_matrix, zh, dual_rhs);
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

  solver.solve(dual_system_matrix, zh, dual_rhs, preconditioner);

  dual_constraints.distribute(zh);

  cout<<"      Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

template <int dim>
void Problem<dim>::output_dual_results() {

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dual_dof_handler);
  data_out.add_data_vector(zh, "zh");
  
  Vector<double> constraint_indicator(dual_dof_handler.n_dofs());
  for (unsigned int i = 0; i < dual_dof_handler.n_dofs(); ++i)
    constraint_indicator(i) = dual_constraints.is_constrained(i) ? 1.0 : 0.0;

  data_out.add_data_vector(constraint_indicator, "constraints");
  data_out.add_data_vector(Rg_dual, "Rg_dual");
  data_out.add_data_vector(uh0_on_dual_space, "uh0_on_dual_space");
  data_out.add_data_vector(Rg_plus_uh0hat, "Rg_plus_uh0hat");
  data_out.add_data_vector(error_indicators, "error_indicators");
  
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
void Problem<dim>::estimate_error(){
	// ------------------------------------------------------------      
  // PROJECTIONS: for both LOCAL and GLOBAL estimates
  // ------------------------------------------------------------      

  uh0_on_dual_space.reinit(dual_dof_handler.n_dofs());
  const auto sol_size = uh0_on_dual_space.size();
  FETools::interpolate(primal_dof_handler,
                        uh0, 
                        dual_dof_handler, 
                        dual_constraints,
                        uh0_on_dual_space);

  Rg_dual.reinit(dual_dof_handler.n_dofs());
  FETools::interpolate(primal_dof_handler,
                        Rg_primal, 
                        dual_dof_handler, 
                        dual_constraints,
                        Rg_dual);

  dual_constraints.distribute(Rg_dual);

  // Subtract from dual solution its projection on the primal solution FE space
  Vector<double> dual_weights(dual_dof_handler.n_dofs());
  FETools::interpolation_difference(dual_dof_handler,
                                    dual_constraints,
                                    zh,
                                    primal_dof_handler,
                                    primal_constraints,
                                    dual_weights);               // zh-∏zh
  dual_constraints.distribute(dual_weights); 

  // ------------------------------------------------------------      
  // RETRIEVE LIFTING: Rg + uh0hat
  // ------------------------------------------------------------      
  // ! Here we were summing up Rg and then forgetting the vector forever

  // uh0^hat + Rg
  Rg_plus_uh0hat.reinit(dual_dof_handler.n_dofs());
  Rg_plus_uh0hat = uh0_on_dual_space;
  Rg_plus_uh0hat += Rg_dual;
  dual_constraints.distribute(Rg_plus_uh0hat);

  // ------------------------------------------------------------      
  // LOCAL ESTIMATE: integrate over cells
  // ------------------------------------------------------------      

  error_indicators.reinit(triangulation.n_active_cells());

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
  // ------------------------------------------------------------      
  
  double dx = 0.0; double lx = 0.0;

  // Compute   u' A z              NOTE: z is actually the dual weights (zh-∏zh)
  Vector<double> temp(dual_dof_handler.n_dofs());
  dual_system_matrix.vmult(temp,dual_weights);    // A z      NB: A is condensed
  dual_constraints.distribute(temp);

  for(unsigned int i=0;i<sol_size;++i)
      dx+=uh0_on_dual_space(i)*temp(i);           // u' A z

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
    fe_values.get_function_gradients(Rg_dual, rg_gradients);

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
  double global_error = std::abs(-lx -dx);


  // ------------------------------------------------------------      
  // SUM UP and OUTPUT
  // ------------------------------------------------------------      

  cout<<"      Global error = " <<  global_error << endl;

  // Sum contribution of each cell's local error to get a global estimate
  double global_error_as_sum_of_cell_errors=0.0;
  for(unsigned int i=0; i<error_indicators.size(); i++)
      global_error_as_sum_of_cell_errors += error_indicators[i];
  global_error_as_sum_of_cell_errors = std::abs(global_error_as_sum_of_cell_errors);
  cout<<"      Global error as sum of cells' errors = " << global_error_as_sum_of_cell_errors << endl;

  // Take absolute value of each error
  for (double &error_indicator : error_indicators) {
      error_indicator = std::fabs(error_indicator);
  }
}

template <int dim>
void Problem<dim>::refine_mesh() {

	GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,error_indicators, 0.6, 0);
  
  // Prepare the solution transfer object
  SolutionTransfer<dim> primal_solution_transfer(primal_dof_handler);
  // take a copy of the solution vector
  Vector<double> old_Rg_dual_dof_values(Rg_primal);
  // Prepare for refinement (older versions of deal.II)
  primal_solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_dual_dof_values);
  // Perform the refinement
  triangulation.execute_coarsening_and_refinement();
  // Reinitialize the DoFHandler for the refined mesh
  primal_dof_handler.distribute_dofs(primal_fe);
  
  
  // Reinitialize Rg_primal to match the new DoF layout after refinement
  Rg_primal.reinit(primal_dof_handler.n_dofs());
  // Transfer the old values to the new DoFs, accounting for hanging nodes
  primal_solution_transfer.interpolate(old_Rg_dual_dof_values, Rg_primal);

  // Handle boundary conditions again (for hanging nodes)
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(primal_dof_handler, 1, exact_solution_function, boundary_values);
  VectorTools::interpolate_boundary_values(primal_dof_handler, 9, exact_solution_function, boundary_values);
  for (const auto &boundary_value : boundary_values)
    Rg_primal(boundary_value.first) = boundary_value.second;
}


template <int dim>
void Problem<dim>::test_convergence(){
  Point<dim> sensor_1(0.00025, 0.0005);
  Point<dim> sensor_2(-0.00025, 0.0005);
  Point<dim> sensor_3(0.0031, 0.0035);

  cout<<"   Convergence test:"<<endl;
  double uh_at_sensor_1 = VectorTools::point_value(primal_dof_handler, uh, sensor_1);
  double uex_at_sensor_1 = exact_solution_function.value(sensor_1);
  std::cout << "      Sensor 1:" << endl
            << "         uh      =  " << uh_at_sensor_1 << endl
            << "         u_ex    =  " << uex_at_sensor_1 << endl
            << "         abs_err =  " << std::fabs(uh_at_sensor_1-uex_at_sensor_1) <<endl;
  
  double uh_at_sensor_2 = VectorTools::point_value(primal_dof_handler, uh, sensor_2);
  double uex_at_sensor_2 = exact_solution_function.value(sensor_2);
  std::cout << "      Sensor 2:" << endl
            << "         uh      =  " << uh_at_sensor_2 << endl
            << "         u_ex    =  " << uex_at_sensor_2 << endl
            << "         abs_err =  " << std::fabs(uh_at_sensor_2-uex_at_sensor_2) <<endl;
  
  double uh_at_sensor_3 = VectorTools::point_value(primal_dof_handler, uh, sensor_3);
  double uex_at_sensor_3 = exact_solution_function.value(sensor_3);
  std::cout << "      Sensor 3:" << endl
            << "         uh      =  " << uh_at_sensor_3 << endl
            << "         u_ex    =  " << uex_at_sensor_3 << endl
            << "         abs_err =  " << std::fabs(uh_at_sensor_3-uex_at_sensor_3) <<endl;
  
 
}


// #######################################
// Template initialization
// #######################################
template class Problem<2>;