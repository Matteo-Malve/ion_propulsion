#include "Problem.h"

template <int dim>
void Problem<dim>::create_mesh()
{
    auto X = constants.X;
    auto L = constants.L;
    auto R = constants.R;
	const std::string filename = "input_mesh.msh";
	cout << "Reading file: " << filename << endl;
	std::ifstream input_file(filename);
	GridIn<2>       grid_in;
	grid_in.attach_triangulation(triangulation);
	grid_in.read_msh(input_file);

	const types::manifold_id emitter = 1;
	EmitterGeometry<2> emitter_manifold(constants);

	for (auto &cell : triangulation.active_cell_iterators()) {
		if (cell->at_boundary()) {
			for (unsigned int f=0; f<4; ++f) {
				if (cell->face(f)->at_boundary()) {

					const Point<dim> fc = cell->face(f)->center();

					if (fc[1] > 0 && fc[1] < 3*R && std::abs(fc[0]-X) < L/2) {
						for (unsigned int v = 0; v < 2; ++v)
							cell->face(f)->vertex(v)[1] =  get_emitter_height(constants.X, constants.L, constants.R, cell->face(f)->vertex(v)[0]);
						cell->face(f)->set_manifold_id(emitter);
						//cout << "Set manifold in " << std::abs(fc[0])-R << endl;
					}
				}
			}
		}
	}

	triangulation.set_all_manifold_ids_on_boundary(1, emitter);
	triangulation.set_manifold(emitter, emitter_manifold);

	cout <<"   Set Manifolds"<< endl;

	for (unsigned int i = 0; i < pre_refinement_steps; ++i) {

		Vector<float> criteria(triangulation.n_active_cells());
		cout  << "Active cells " << triangulation.n_active_cells() << endl;
		unsigned int ctr = 0;

		for (auto &cell : triangulation.active_cell_iterators()) {

			const Point<dim> c = cell->center();
			const double d = std::sqrt( (c[0]-X)*(c[0]-X) + c[1]*c[1]);

			if ( d <= pre_refinement_steps/(i+1)*2.*R)
				criteria[ctr++] = 1;
			else
				criteria[ctr++] = 0;
		}
		GridRefinement::refine(triangulation, criteria, 0.5);
		triangulation.execute_coarsening_and_refinement();
	}

	// Refine twice near the edges: only for custom meshes, if needed
	/*Vector<float> criteria(triangulation.n_active_cells());
	cout  << "Active cells " << triangulation.n_active_cells() << endl;
	unsigned int ctr = 0;

	const double left_edge = X - L/2.;
	const double right_edge = X + L/2.;

	for (auto &cell : triangulation.active_cell_iterators()) {

		const Point<dim> c = cell->center();

		if ( c[1] < 2.5e-5 && ( (abs(c[0]-left_edge) <= 1.e-5) ||  (abs(c[0]-right_edge) <= 1.e-5) ) )
			criteria[ctr++] = 1;
		else
			criteria[ctr++] = 0;
	}
	GridRefinement::refine(triangulation, criteria, 0.5);
	triangulation.execute_coarsening_and_refinement();*/

	cout  << "Final number of active cells: " << triangulation.n_active_cells() << endl;
	cout <<"   Executed Pre-Refinement"<< endl;

}


template <int dim>
void Problem<dim>::setup_primal_system()
{
    primal_dof_handler.distribute_dofs(primal_fe);

    primal_constraints.clear();
	DoFTools::make_hanging_node_constraints(primal_dof_handler, primal_constraints);
	primal_constraints.close();

    DynamicSparsityPattern dsp(primal_dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(primal_dof_handler, dsp, primal_constraints, false);
    primal_sparsity_pattern.copy_from(dsp);

    primal_system_matrix.reinit(primal_sparsity_pattern);

    uh0.reinit(primal_dof_handler.n_dofs());
    primal_solution.reinit(primal_dof_handler.n_dofs());
    primal_rhs.reinit(primal_dof_handler.n_dofs());

    cout << "Primal problem DoFs: " << primal_dof_handler.n_dofs() << endl;
}

template <int dim>
void Problem<dim>::setup_dual_system()
{
    dual_dof_handler.distribute_dofs(dual_fe);

    dual_constraints.clear();
	DoFTools::make_hanging_node_constraints(dual_dof_handler, dual_constraints);
	dual_constraints.close();

    DynamicSparsityPattern dsp(dual_dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dual_dof_handler, dsp, dual_constraints, false);
    dual_sparsity_pattern.copy_from(dsp);

    dual_system_matrix.reinit(dual_sparsity_pattern);

    dual_solution.reinit(dual_dof_handler.n_dofs());
    dual_rhs.reinit(dual_dof_handler.n_dofs());

    cout << "Dual problem DoFs: " << dual_dof_handler.n_dofs() << endl;
}


template <int dim>
void Problem<dim>::assemble_primal_system()
{
    auto eps_0 = constants.eps_0;
    auto eps_r = constants.eps_r;
	const QGauss <dim> quadrature(dual_fe.degree + 1);
    FEValues<dim> fe_values(primal_fe,
                            quadrature,
                            update_values | update_gradients | update_quadrature_points |
                            update_JxW_values);

    const unsigned int dofs_per_cell = primal_fe.n_dofs_per_cell();
    const unsigned int n_q_points = quadrature.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double> rhs_values(n_q_points);

    QGauss<dim>          Rg_quadrature(2);
    FEValues<dim>        Rg_fe_values(primal_fe,
                                   Rg_quadrature,
                                   update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values);
    const unsigned int Rg_n_q_points = Rg_quadrature.size();

    for (const auto &cell : primal_dof_handler.active_cell_iterators())
    {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

        // Compute A_loc
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
            for (const unsigned int i : fe_values.dof_indices())
                for (const unsigned int j : fe_values.dof_indices())
                    cell_matrix(i, j) += eps_r*eps_0*
                            (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                             fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                             fe_values.JxW(q_index));           // dx
        }

        // Compute RHS
		rhs_function.value_list(fe_values.get_quadrature_points(),rhs_values);
		// Compute f_loc
		for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
			for (unsigned int i = 0; i < dofs_per_cell; ++i)
				cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
								rhs_values[q_point] *               // f(x_q)
								fe_values.JxW(q_point));            // dx

		Rg_fe_values.reinit(cell);
		// Cycle over quadrature nodes
		for (unsigned int q_point = 0; q_point < Rg_n_q_points; ++q_point){
			// evaluate grad_Rg
			auto & quad_point_coords = Rg_fe_values.get_quadrature_points()[q_point];
			auto grad_Rg_xq = evaluate_grad_Rg(constants.X, constants.L, constants.R, constants.Vmax, quad_point_coords[0],quad_point_coords[1]);

			// assemble A_loc(Rg,v)
			for (unsigned int i = 0; i < dofs_per_cell; ++i)
				cell_rhs(i) -= eps_r*eps_0*
								(Rg_fe_values.shape_grad(i, q_point) *     // grad phi_i(x_q)
								grad_Rg_xq *                               // grad_Rg(x_q)
								Rg_fe_values.JxW(q_point));                // dx

		}

        // Local to global
        cell->get_dof_indices(local_dof_indices);
        primal_constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, primal_system_matrix, primal_rhs);
    }

    // Apply boundary values
    {
      std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

      VectorTools::interpolate_boundary_values(primal_dof_handler,1, Functions::ConstantFunction<dim>(0.), emitter_boundary_values);
      MatrixTools::apply_boundary_values(emitter_boundary_values, primal_system_matrix, primal_solution, primal_rhs);

      VectorTools::interpolate_boundary_values(primal_dof_handler,2, Functions::ConstantFunction<dim>(0.), collector_boundary_values);
      MatrixTools::apply_boundary_values(collector_boundary_values, primal_system_matrix, primal_solution, primal_rhs);
    }

    cout<<"   Assembled primal system" <<endl;
}

template <int dim>
void Problem<dim>::assemble_dual_system()
{
    auto eps_0 = constants.eps_0;
    auto eps_r = constants.eps_r;
    FEValues<dim> fe_values(dual_fe,
                            dual_quadrature,
                            update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = dual_fe.n_dofs_per_cell();
    const unsigned int n_q_points = dual_quadrature.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dual_dof_handler.active_cell_iterators())
    {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

        // Compute A_loc
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
            for (const unsigned int i : fe_values.dof_indices())
                for (const unsigned int j : fe_values.dof_indices())
                    cell_matrix(i, j) += eps_r*eps_0*
                            (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                             fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                             fe_values.JxW(q_index));           // dx
        }

        // Compute RHS
		for (const auto &face: cell->face_iterators()) {
			if (face->at_boundary()) {
				if (face->boundary_id() == 1) {
					cell_rhs = 0;
					fe_values.reinit(cell);

					// Retrieve the components of the normal vector to the boundary face
					Tensor<1, dim> n = emitter_normal(constants.X, constants.L, constants.R, face->center());

					// Compute the flux for this cell
					for (unsigned int q = 0; q < n_q_points; ++q) {
						for (unsigned int i = 0; i < dofs_per_cell; ++i) {
							for (unsigned int k = 0; k < dim; ++k) {
								cell_rhs[i] += fe_values.shape_grad(i, q)[k] * (-n[k]);
							}
							cell_rhs[i] *= fe_values.JxW(q);
						}
					}
				}
			}
		}

        // Local to global
        cell->get_dof_indices(local_dof_indices);
        dual_constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, dual_system_matrix, dual_rhs);
    }

    // Apply boundary values
    {
      std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

      VectorTools::interpolate_boundary_values(dual_dof_handler,1, Functions::ConstantFunction<dim>(0.), emitter_boundary_values);
      MatrixTools::apply_boundary_values(emitter_boundary_values, dual_system_matrix, dual_solution, dual_rhs);

      VectorTools::interpolate_boundary_values(dual_dof_handler,2, Functions::ConstantFunction<dim>(0.), collector_boundary_values);
      MatrixTools::apply_boundary_values(collector_boundary_values, dual_system_matrix, dual_solution, dual_rhs);
    }

    cout<<"   Assembled dual system" <<endl;
}


template <int dim>
void Problem<dim>::solve_primal()
{
  const unsigned int it_max = 1e+4;
  const double rel_tol = 1.e-6*primal_rhs.l2_norm();
  const double abs_tol = 1.e-12;

  const double tol = abs_tol + rel_tol;
  SolverControl            solver_control(it_max, tol);
  SolverCG<Vector<double>> solver(solver_control);

  double relaxation_parameter = 1.2;
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(primal_system_matrix, relaxation_parameter);

  solver.solve(primal_system_matrix, primal_solution, primal_rhs, preconditioner);

  primal_constraints.distribute(primal_solution);

  cout<<"   Solved primal problem: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

template <int dim>
void Problem<dim>::solve_dual()
{
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

  cout<<"   Solved dual problem: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}


template <int dim>
void Problem<dim>::output_primal_results()
{
    auto g = constants.g;
    auto nn = constants.nn;
    const Point<dim> evaluation_point(0.5*g, 0.1*g);
    const double x = VectorTools::point_value(primal_dof_handler, primal_solution, evaluation_point);
    cout << "   Potential at sample point (" << evaluation_point[0] << "," << evaluation_point[1] << "): " << x << endl;

	Gradient el_field(constants);
	IonizationArea ion_area(constants);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(primal_dof_handler);
    data_out.add_data_vector(primal_solution, "Potential");
    data_out.add_data_vector(primal_solution, el_field);
	data_out.add_data_vector(primal_solution, ion_area);
    data_out.build_patches(); // mapping

    std::string filename;
	filename = Utilities::int_to_string(nn, 1) + "R_test_solution-" + Utilities::int_to_string(cycle, 1) + ".vtk";
    DataOutBase::VtkFlags vtk_flags;
    vtk_flags.compression_level = DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
    data_out.set_flags(vtk_flags);
    std::ofstream output(filename);
    data_out.write_vtk(output);
}

template <int dim>
void Problem<dim>::output_dual_results()
{
    auto nn = constants.nn;
	Gradient el_field(constants);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dual_dof_handler);
    data_out.add_data_vector(dual_solution, "Potential");
    data_out.add_data_vector(dual_solution, el_field);
    data_out.build_patches(); // mapping

    std::string filename;
	filename =  Utilities::int_to_string(nn, 1) + "R_dual_test_solution-" + Utilities::int_to_string(cycle, 1) + ".vtk";
    DataOutBase::VtkFlags vtk_flags;
    vtk_flags.compression_level = DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
    data_out.set_flags(vtk_flags);
    std::ofstream output(filename);
    data_out.write_vtk(output);
}


template <int dim>
void Problem<dim>::refine_mesh()
{
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
    cout<<"   Global error = " <<  global_error << endl
        <<"   Global error as sum of cells' errors = " << global_error_as_sum_of_cell_errors << endl;

    // Take absolute value of each error
    for (float &error_indicator : error_indicators) {
        error_indicator = std::fabs(error_indicator);
    }

	GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,error_indicators, 0.8, 0.0); // CHANGED FROM: 0.8, 0.02

    // Execute refinement
    triangulation.execute_coarsening_and_refinement();
    cout<<"   Executed coarsening and refinement"<<endl;
}

template <int dim>
double Problem<dim>::estimate_error(Vector<float> &error_indicators) const
{
    auto eps_0 = constants.eps_0;
    auto eps_r = constants.eps_r;
    Vector<double> interpolated_primal_solution(dual_dof_handler.n_dofs());
    FETools::interpolate(primal_dof_handler,uh0, dual_dof_handler, dual_constraints,interpolated_primal_solution);

    dual_constraints.distribute(interpolated_primal_solution); // ADDED

    // Subtract from dual solution its projection on the primal solution FE space
    Vector<double> dual_weights(dual_dof_handler.n_dofs());
    FETools::interpolation_difference(dual_dof_handler,
                                      dual_constraints,
                                      dual_solution,
                                      primal_dof_handler,
                                      primal_constraints,
                                      dual_weights);                            // zh-∏zh

    const auto sol_size = interpolated_primal_solution.size();

    // 1) COMPUTE ERROR INDICATORS
    // Compute Rg_plus_uh0hat
    Vector<double> Rg_plus_uh0hat(sol_size);
    Vector<double> Rg_vector(sol_size);
    VectorTools::interpolate(dual_dof_handler,Evaluate_Rg<dim>(constants),Rg_vector);

    dual_constraints.distribute(Rg_vector); // ADDED

    for(unsigned int i=0;i<sol_size;i++)
        Rg_plus_uh0hat(i) = Rg_vector(i) + interpolated_primal_solution(i);

    dual_constraints.distribute(Rg_plus_uh0hat); // ADDED

    // INTEGRATE OVER CELLS
    FEValues<dim> fe_values(dual_fe,
                            dual_quadrature,
                            update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned int n_q_points = dual_quadrature.size();

    std::vector<Tensor<1,dim>> cell_primal_gradients(n_q_points);
    std::vector<Tensor<1,dim>> cell_dual_gradients(n_q_points);

    double sum;

    for (const auto &cell : dual_dof_handler.active_cell_iterators()){

        fe_values.reinit(cell);
        fe_values.get_function_gradients(interpolated_primal_solution, cell_primal_gradients);
        fe_values.get_function_gradients(dual_weights, cell_dual_gradients);

        // Numerically approximate the integral of the scalar product between the gradients of the two
        sum = 0;
        for (unsigned int p = 0; p < n_q_points; ++p) {
            sum +=  eps_r*eps_0*
                    ((cell_primal_gradients[p] * cell_dual_gradients[p]  )   // Scalar product btw Tensors
                     * fe_values.JxW(p));
        }
        error_indicators(cell->active_cell_index()) += (0. - sum);

    }


    // 2) EVALUATE GLOBAL ERROR
    double dx = 0.0; double lx = 0.0;

    // Compute   u' A z              ;NOTE: z is actually the dual weights (zh-∏zh)
    Vector<double> temp(dual_dof_handler.n_dofs());
    dual_system_matrix.vmult(temp,dual_weights);    // A z

    for(unsigned int i=0;i<sol_size;++i)
        dx+=interpolated_primal_solution(i)*temp(i);           // u' A z

    // Compute   a(Rg,φ)
    const unsigned int dofs_per_cell = dual_fe.n_dofs_per_cell();

    Vector<double> F(dual_dof_handler.n_dofs());
    Vector<double> cell_F(dofs_per_cell);
    std::vector <types::global_dof_index> local_dof_indices(dofs_per_cell);

    QGauss<dim>          Rg_quadrature(5);
    FEValues<dim>        Rg_fe_values(dual_fe,
                                      Rg_quadrature,
                                      update_values | update_gradients | update_quadrature_points |
                                      update_JxW_values);
    const unsigned int Rg_n_q_points = Rg_quadrature.size();

    for (const auto &cell: dual_dof_handler.active_cell_iterators()) {
        cell_F = 0;
        Rg_fe_values.reinit(cell);
        for (unsigned int q_point = 0; q_point < Rg_n_q_points; ++q_point){
            // evaluate_grad_Rg
            auto & quad_point_coords = Rg_fe_values.get_quadrature_points()[q_point];
            auto grad_Rg_xq = evaluate_grad_Rg(constants.X, constants.L, constants.R, constants.Vmax, quad_point_coords[0],quad_point_coords[1]);
            // assemble a(Rg,φ)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_F(i) += eps_r*eps_0*
                               (Rg_fe_values.shape_grad(i, q_point) *      // grad phi_i(x_q)
                                grad_Rg_xq *                               // grad_Rg(x_q)
                                Rg_fe_values.JxW(q_point));                // dx
        }
        // Local to Global
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            F(local_dof_indices[i]) += cell_F(i);

        //dual_constraints.distribute_local_to_global(cell_F,local_dof_indices,F);
    }

    // Compute     z' F    or    z•F
    for(unsigned int i=0;i<dual_weights.size();i++)
        lx += dual_weights(i) * F(i);

    // Return             η = |r(z)|  = | - z' F - u' A z |
    // or, precisely,     η = |r(zk-∏zk)|  = | - (zk-∏zk)' F - uh0' A (zk-∏zk) |
    return std::abs(-lx -dx);

}


template <int dim>
void Problem<dim>::run()
{
	create_mesh();

	while (cycle <= max_refinements) {

		setup_primal_system();  // make
		assemble_primal_system();       // condense
		solve_primal();         // distribute
		uh0 = primal_solution;

		// Retrieve lifting
		{
			Vector<double> Rg_vector(primal_solution.size());
			VectorTools::interpolate(primal_dof_handler,Evaluate_Rg<dim>(constants),Rg_vector);
			primal_constraints.distribute(Rg_vector);       // distribute
			primal_solution += Rg_vector;               // uh = u0 + Rg
			primal_constraints.distribute(primal_solution);     // distribute
		}
		output_primal_results();

		setup_dual_system();
		assemble_dual_system();
		solve_dual();
		output_dual_results();

		refine_mesh();

		cout << " 	Elapsed CPU time: " << timer.cpu_time() << " seconds for " + Utilities::int_to_string(cycle,1) <<" cycles" << endl;

		++cycle;
	}
}



// #######################################
// Template initialization
// #######################################
template class Problem<2>;

