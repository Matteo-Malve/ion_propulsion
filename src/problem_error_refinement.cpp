#include "problem.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

using namespace dealii;
using std::cout;
using std::endl;
namespace plt = matplotlibcpp;

// -----------------------------------------
// GOAL ORIENTED REFINEMENT
// -----------------------------------------

template <int dim>
void Problem<dim>::interpolate_between_primal_and_dual(){
  // ------------------------------------------------------------      
  // PROJECTIONS: for both LOCAL and GLOBAL estimates
  // ------------------------------------------------------------      

  uh0_on_dual_space.reinit(dual_dof_handler.n_dofs());
  FETools::interpolate(primal_dof_handler,
                        uh0, 
                        dual_dof_handler, 
                        dual_constraints,
                        uh0_on_dual_space);
  dual_constraints.distribute(uh0_on_dual_space); // just added, shouldn't have been necessary

  Rg_dual.reinit(dual_dof_handler.n_dofs());
  FETools::interpolate(primal_dof_handler,
                        Rg_primal, 
                        dual_dof_handler, 
                        dual_constraints,
                        Rg_dual);

  dual_constraints.distribute(Rg_dual);

  // Subtract from dual solution its projection on the primal solution FE space
  dual_weights.reinit(dual_dof_handler.n_dofs());
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

  Rg_plus_uh0hat.reinit(dual_dof_handler.n_dofs());
  Rg_plus_uh0hat = uh0_on_dual_space;
  Rg_plus_uh0hat += Rg_dual;
  dual_constraints.distribute(Rg_plus_uh0hat);
}

template <int dim>
void Problem<dim>::local_estimate(){
  // ------------------------------------------------------------      
  // LOCAL ESTIMATE: integrate over cells
  // ------------------------------------------------------------      

  error_indicators.reinit(triangulation.n_active_cells());

  const QGauss<dim> quadrature(dual_dof_handler.get_fe().degree + 1);
  FEValues<dim> fe_values(dual_fe,
                          quadrature,
                          update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = quadrature.size();

  std::vector<Tensor<1,dim>> cell_Rg_plus_uh0hat_gradients(n_q_points);
  std::vector<Tensor<1,dim>> cell_dual_weights_gradients(n_q_points);
  std::vector<double> cell_rhs_values(n_q_points);
  std::vector<double> cell_dual_weights(n_q_points);

  double sum;

  for (const auto &cell : dual_dof_handler.active_cell_iterators()){

    fe_values.reinit(cell);
    fe_values.get_function_gradients(Rg_plus_uh0hat, cell_Rg_plus_uh0hat_gradients);
    fe_values.get_function_gradients(dual_weights, cell_dual_weights_gradients);
    fe_values.get_function_values(dual_weights, cell_dual_weights);
    rhs_function.value_list(fe_values.get_quadrature_points(), cell_rhs_values);

    // Numerically approximate the integral of the scalar product between the gradients of the two
    sum = 0;
    
    for (unsigned int p = 0; p < n_q_points; ++p) {
        sum +=  eps_r * eps_0 *                        
                ((cell_Rg_plus_uh0hat_gradients[p] * cell_dual_weights_gradients[p]  )   // Scalar product btw Tensors
                  * fe_values.JxW(p));
        sum += (cell_rhs_values[p] * cell_dual_weights[p] * fe_values.JxW(p));
    }
    error_indicators(cell->active_cell_index()) -= sum;  

  }
}

template <int dim>
void Problem<dim>::local_estimate_face_jumps(){
  error_indicators_face_jumps.reinit(triangulation.n_active_cells());
  // Having computed the dual weights we now proceed with computing the cell and face residuals 
  // of the primal solution. First we set up a map between face iterators and their jump term 
  // contributions of faces to the error estimator. The reason is that we compute the jump terms 
  // only once, from one side of the face, and want to collect them only afterwards when looping 
  // over all cells a second time.
  // We initialize this map already with a value of -1e20 for all faces, since this value will 
  // stand out in the results if something should go wrong and we fail to compute the value for a 
  // face for some reason. Secondly, this initialization already makes the std::map object allocate 
  // all objects it may possibly need.

  // Map between face iterators and their jump term contributions
  std::map<typename DoFHandler<dim>::face_iterator, double> face_integrals;
  
  // Quadrature rule for face integration
  const QGauss<dim> quadrature(dual_fe.degree + 1);
  const QGauss<dim - 1> face_quadrature(dual_fe.degree + 1);

  const unsigned int n_q_points = quadrature.size();
  const unsigned int face_n_q_points = face_quadrature.size();

  std::vector<double> cell_rhs_values(n_q_points);
  std::vector<double> cell_dual_weights(n_q_points);
  std::vector<double> cell_laplacians(n_q_points);

  // FEValues objects for face integration
  FEValues<dim> fe_values (dual_fe, quadrature,
                                      update_values | update_hessians | 
                                      update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values_current_cell(dual_fe, face_quadrature,
                                   update_values | update_gradients |
                                   update_normal_vectors | update_JxW_values);
  FEFaceValues<dim> fe_face_values_neighbor(dual_fe, face_quadrature,
                                            update_values | update_gradients | 
                                            update_normal_vectors | update_JxW_values);
  FESubfaceValues<dim> fe_subface_values_current_cell(dual_fe, face_quadrature,
                                         update_values | update_gradients | 
                                         update_normal_vectors | update_JxW_values);
  
  // Vectors to store gradients on the face quadrature points
  std::vector<Tensor<1, dim>> cell_gradients(face_n_q_points);
  std::vector<Tensor<1, dim>> neighbor_gradients(face_n_q_points);
  std::vector<double> face_dual_weights(face_n_q_points);

  // Initializing face_integrals to detect any uncomputed values
  for (const auto &cell : dual_dof_handler.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      face_integrals[face] = -1e20;


  for (const auto &cell : dual_dof_handler.active_cell_iterators()){
    // INTEGRATE OVER CELL (with laplacian)
    fe_values.reinit(cell);
    rhs_function.value_list(fe_values.get_quadrature_points(), cell_rhs_values);
    fe_values.get_function_laplacians(Rg_plus_uh0hat, cell_laplacians);
    fe_values.get_function_values(dual_weights, cell_dual_weights);
    double sum = 0;
    
    for (unsigned int p = 0; p < n_q_points; ++p)
      sum += ((cell_rhs_values[p] + cell_laplacians[p]) *
             cell_dual_weights[p] * fe_values.JxW(p));

    error_indicators_face_jumps(cell->active_cell_index()) += sum;
  }

  // LOOP { ESTIMATE ON ONE CELL }
  for (const auto &cell : dual_dof_handler.active_cell_iterators()){
    for (const unsigned int face_no : cell->face_indices()){
      // If this face is part of the boundary, then there is nothing to do
      if (cell->face(face_no)->at_boundary()){
        face_integrals[cell->face(face_no)] = 0;
        continue; 
      }
      // First, if the neighboring cell is on the same level as this one, i.e. neither further 
      // refined not coarser, then the one with the lower index within this level does the work. 
      // In other words: if the other one has a lower index, then skip work on this face:
      if ((cell->neighbor(face_no)->has_children() == false) &&
          (cell->neighbor(face_no)->level() == cell->level()) &&
          (cell->neighbor(face_no)->index() < cell->index()))
        continue;
      // Likewise, we always work from the coarser cell if this and its neighbor differ in refinement. 
      // Thus, if the neighboring cell is less refined than the present one, then do
      // nothing since we integrate over the subfaces when we visit the coarse cell.
      if (cell->at_boundary(face_no) == false)
        if (cell->neighbor(face_no)->level() < cell->level())
          continue;

      // Now we know that we are in charge here, so actually compute the face jump terms.
      if (cell->face(face_no)->has_children() == false){
        // INTEGRATE OVER REGULAR FACE
        std::vector<double> jump_residual(face_n_q_points);

        fe_face_values_current_cell.reinit(cell, face_no);
        fe_face_values_current_cell.get_function_gradients(Rg_plus_uh0hat, cell_gradients);

        fe_face_values_neighbor.reinit(cell->neighbor(face_no), cell->neighbor_of_neighbor(face_no));
        fe_face_values_neighbor.get_function_gradients(Rg_plus_uh0hat, neighbor_gradients);

        Assert(cell->neighbor(face_no).state() == IteratorState::valid,ExcInternalError());

        for (unsigned int q = 0; q < face_n_q_points; ++q)
          jump_residual[q] = ((cell_gradients[q] - neighbor_gradients[q]) 
                              * fe_face_values_current_cell.normal_vector(q));

        fe_face_values_current_cell.get_function_values(dual_weights,face_dual_weights);

        double face_integral = 0;
        for (unsigned int q = 0; q < face_n_q_points; ++q)
          face_integral += (jump_residual[q] * face_dual_weights[q] 
                            * fe_face_values_current_cell.JxW(q));

        Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(),ExcInternalError());
        Assert(face_integrals[cell->face(face_no)] == -1e20, ExcInternalError());
        
        face_integrals[cell->face(face_no)] = face_integral;

      } else {

        // INTEGRATE OVER IRREGULAR FACE

        std::vector<double> jump_residual(face_n_q_points);

        const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
        const unsigned int neighbor_neighbor = cell->neighbor_of_neighbor(face_no);
        const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);

        Assert(neighbor.state() == IteratorState::valid, ExcInternalError());
        Assert(neighbor->has_children(), ExcInternalError());

        for (unsigned int subface_no = 0; subface_no < face->n_children(); ++subface_no){
          const typename DoFHandler<dim>::active_cell_iterator  neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);
          Assert(neighbor_child->face(neighbor_neighbor) == cell->face(face_no)->child(subface_no), ExcInternalError());
          
          fe_subface_values_current_cell.reinit(cell, face_no, subface_no);
          fe_subface_values_current_cell.get_function_gradients(Rg_plus_uh0hat, cell_gradients);

          fe_face_values_neighbor.reinit(neighbor_child, neighbor_neighbor);
          fe_face_values_neighbor.get_function_gradients(Rg_plus_uh0hat, neighbor_gradients);

          for (unsigned int p = 0; p < face_n_q_points; ++p)
            jump_residual[p] = ((neighbor_gradients[p] - cell_gradients[p]) *
                                fe_face_values_neighbor.normal_vector(p));

          fe_face_values_neighbor.get_function_values(dual_weights, face_dual_weights);           

          double face_integral = 0;
          for (unsigned int p = 0; p < face_n_q_points; ++p)
            face_integral += (jump_residual[p] * face_dual_weights[p] *
                              fe_face_values_neighbor.JxW(p));
          
          face_integrals[neighbor_child->face(neighbor_neighbor)] = face_integral;           

        }

        double sum = 0;
        for (unsigned int subface_no = 0; subface_no < face->n_children(); ++subface_no){
          Assert(face_integrals.find(face->child(subface_no)) != face_integrals.end(), ExcInternalError());
          Assert(face_integrals[face->child(subface_no)] != -1e20, ExcInternalError());
          sum += face_integrals[face->child(subface_no)];
        }
        face_integrals[face] = sum;
      }
    }
  }

  unsigned int present_cell = 0;
  for (const auto &cell : dual_dof_handler.active_cell_iterators()){
    for (const auto &face : cell->face_iterators()){
      Assert(face_integrals.find(face) != face_integrals.end(), ExcInternalError());
      error_indicators_face_jumps(present_cell) -= 0.5 * face_integrals[face];
    }
    ++present_cell;
  }
  /*std::cout << "   Estimated error [face jumps]: "
            << std::accumulate(error_indicators_face_jumps.begin(),error_indicators_face_jumps.end(),0.)
            << std::endl;*/

}

template <int dim>
double Problem<dim>::global_estimate(){
  // ------------------------------------------------------------      
  // GLOBAL ESTIMATE
  // ------------------------------------------------------------      
  const QGauss<dim> quadrature(dual_dof_handler.get_fe().degree + 1);
  FEValues<dim> fe_values(dual_fe,
                          quadrature,
                          update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = quadrature.size();
  const auto sol_size = uh0_on_dual_space.size();


  dual_dof_handler.distribute_dofs(dual_fe);
  DynamicSparsityPattern dsp(dual_dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dual_dof_handler, dsp);
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);

  // Globals
  SparseMatrix<double> dual_system_matrix_without_constraints;
  dual_system_matrix_without_constraints.reinit(sparsity_pattern);
  Vector<double> F(dual_dof_handler.n_dofs());

  // FE
  const unsigned int dofs_per_cell = dual_fe.n_dofs_per_cell();

  std::vector<Tensor<1, dim>> rg_gradients(n_q_points);

  // Locals
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_F(dofs_per_cell);
  std::vector<double> cell_rhs_values(quadrature.size());

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  

  for (const auto &cell : dual_dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_F = 0;

    fe_values.get_function_gradients(Rg_dual, rg_gradients);
    rhs_function.value_list(fe_values.get_quadrature_points(), cell_rhs_values);
    
    // Compute A_loc
    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      for (const unsigned int i : fe_values.dof_indices()){
        cell_F(i) += eps_r * eps_0*
                        (fe_values.shape_grad(i, q_index) *     // grad phi_i(x_q)
                        rg_gradients[q_index] *                 // grad_Rg(x_q)
                        fe_values.JxW(q_index));                // dx
        cell_F(i) += (cell_rhs_values[q_index] * fe_values.shape_value(i, q_index) * fe_values.JxW(q_index));                
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) += eps_r*eps_0*
                  (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                    fe_values.JxW(q_index));           // dx
      }
    // Local to global
    cell->get_dof_indices(local_dof_indices);
    for (const unsigned int i : fe_values.dof_indices()){
      F(local_dof_indices[i]) += cell_F(i);
      for (const unsigned int j : fe_values.dof_indices())
        dual_system_matrix_without_constraints.add(local_dof_indices[i],        // <-- dual_system_matrix_without_constraints
                                                    local_dof_indices[j],
                                                    cell_matrix(i, j)); 
    }
  }

  double dx = 0.0; double lx = 0.0;

  Vector<double> temp(dual_dof_handler.n_dofs());
  dual_system_matrix_without_constraints.vmult(temp,dual_weights);    // A z 

  for(unsigned int i=0;i<sol_size;++i)
      dx+=uh0_on_dual_space(i)*temp(i);           // u' A z   // Prova anche a togliere distribute da uh0_on_dual_space

  // Compute     z' F    or    z•F
  for(unsigned int i=0;i<dual_weights.size();i++)
    lx += dual_weights(i) * F(i);

  // Return             η = |r(z)|  = | - z' F - u' A z |
  // or, precisely,     η = |r(zk-∏zk)|  = | - (zk-∏zk)' F - uh0' A (zk-∏zk) |
  double global_error = std::abs(-lx -dx);
  return global_error;
}


template <int dim>
void Problem<dim>::estimate_error(){
  // ------------------------------------------------------------      
  // CALL all steps
  // ------------------------------------------------------------   
  
  interpolate_between_primal_and_dual();
  local_estimate();
  local_estimate_face_jumps();
  double global_error = global_estimate();

  // ------------------------------------------------------------      
  // SUM UP and OUTPUT
  // ------------------------------------------------------------      
  
  cout<<"      Global error = " <<  global_error << endl;
  goal_oriented_global_errors.push_back(global_error);

  // Sum contribution of each cell's local error to get a global estimate
  double global_error_as_sum_of_cell_errors=0.0;
  for(unsigned int i=0; i<error_indicators.size(); i++)
      global_error_as_sum_of_cell_errors += error_indicators[i];
  global_error_as_sum_of_cell_errors = std::abs(global_error_as_sum_of_cell_errors);
  cout<<"      Global error as sum of cells' errors = " << global_error_as_sum_of_cell_errors << endl;
  goal_oriented_local_errors.push_back(global_error_as_sum_of_cell_errors);

  // Sum contribution of each cell's local error to get a global estimate
  double global_error_as_sum_of_cell_errors_face_jumps = 0.0;
  for(unsigned int i=0; i<error_indicators_face_jumps.size(); i++)
      global_error_as_sum_of_cell_errors_face_jumps += error_indicators_face_jumps[i];
  global_error_as_sum_of_cell_errors_face_jumps = std::abs(global_error_as_sum_of_cell_errors_face_jumps);
  cout<<"      Global error as sum of cells' errors [face jumps] = " << global_error_as_sum_of_cell_errors_face_jumps << endl;
  goal_oriented_local_errors_face_jumps.push_back(global_error_as_sum_of_cell_errors_face_jumps);

  // Take absolute value of each error
  for (double &error_indicator : error_indicators)
    error_indicator = std::fabs(error_indicator);
  
  for (double &error_indicator : error_indicators_face_jumps) 
    error_indicator = std::fabs(error_indicator);

  // ------------------------------------------------------------      
  // PLOT
  // ------------------------------------------------------------      

  plt::figure_size(800, 600);
  plt::clf(); // Clear previous plot

  plt::named_semilogy("local_estimate", cycles, goal_oriented_local_errors, "b-o");
  plt::named_semilogy("local_estimate_face_jumps", cycles, goal_oriented_local_errors_face_jumps, "b--o");
  plt::named_semilogy("global_estimate", cycles, goal_oriented_global_errors, "r-o");

  plt::xlabel("Cycle");
  plt::ylabel("Error");
  plt::title("Goal-oriented error");
  plt::legend();
  plt::grid(true);

  // Save the plot
  plt::save(TEST_NAME+"-goal_oriented_convergence.png");

  // ------------------------------------------------------------      
  // TABLE
  // ------------------------------------------------------------      

  GO_table.add_value("cycle", cycle);
  GO_table.add_value("cells", triangulation.n_active_cells());
  GO_table.add_value("global", global_error);
  GO_table.add_value("local", global_error_as_sum_of_cell_errors);
  GO_table.add_value("l. jumps", global_error_as_sum_of_cell_errors_face_jumps);

}

template <int dim>
void Problem<dim>::refine_mesh() {
  if(REFINEMENT_STRATEGY == "GO"){
    //GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,error_indicators, 0.6, 0.01);
    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,error_indicators_face_jumps, 0.6, 0.1);
    //GridRefinement::refine_and_coarsen_fixed_number(triangulation,error_indicators_face_jumps, 0.05, 0);
    triangulation.prepare_coarsening_and_refinement();
    // Prepare the solution transfer object
    SolutionTransfer<dim> primal_solution_transfer(primal_dof_handler);
    // take a copy of the solution vector
    Vector<double> old_Rg_dof_values(Rg_primal);
    // Prepare for refinement (older versions of deal.II)
    primal_solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_dof_values);
    // Perform the refinement
    triangulation.execute_coarsening_and_refinement();
    // Reinitialize the DoFHandler for the refined mesh
    primal_dof_handler.distribute_dofs(primal_fe);
    
    
    // Reinitialize Rg_primal to match the new DoF layout after refinement
    Rg_primal.reinit(primal_dof_handler.n_dofs());
    // Transfer the old values to the new DoFs, accounting for hanging nodes
    primal_solution_transfer.interpolate(old_Rg_dof_values, Rg_primal);

    // Handle boundary conditions again (for hanging nodes)
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(primal_dof_handler, 1, exact_solution_function, boundary_values);
    VectorTools::interpolate_boundary_values(primal_dof_handler, 9, exact_solution_function, boundary_values);
    for (const auto &boundary_value : boundary_values)
      Rg_primal(boundary_value.first) = boundary_value.second;

  } else if(REFINEMENT_STRATEGY == "GlobRef"){

    Vector<double> old_Rg_values = Rg_primal;

    SolutionTransfer<dim> solution_transfer(primal_dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(old_Rg_values);

    triangulation.refine_global(1);

    primal_dof_handler.distribute_dofs(primal_fe);

    Rg_primal.reinit(primal_dof_handler.n_dofs());
    solution_transfer.interpolate(old_Rg_values, Rg_primal);

  } else {
    cout << "Refinement strategy undefined" << endl;
    abort();
  }
}


// #######################################
// Template initialization
// #######################################
template class Problem<2>;