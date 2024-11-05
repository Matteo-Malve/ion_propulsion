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
}

template <int dim>
void Problem<dim>::local_estimate_face_jumps(){
 
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
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  

  for (const auto &cell : dual_dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_F = 0;

    fe_values.get_function_gradients(Rg_dual, rg_gradients);

    // Compute A_loc
    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      for (const unsigned int i : fe_values.dof_indices()){
        cell_F(i) += eps_r * eps_0*
                        (fe_values.shape_grad(i, q_index) *     // grad phi_i(x_q)
                        rg_gradients[q_index] *                 // grad_Rg(x_q)
                        fe_values.JxW(q_index));                // dx
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
  //local_estimate_face_jumps();
  double global_error = global_estimate();

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
  goal_oriented_local_errors.push_back(global_error_as_sum_of_cell_errors);

  // Take absolute value of each error
  for (double &error_indicator : error_indicators) {
      error_indicator = std::fabs(error_indicator);
  }

  // Plot data
  plt::figure_size(800, 600);
  plt::clf(); // Clear previous plot

  plt::named_semilogy("local_estimate", cycles, goal_oriented_local_errors, "b-o");

  plt::xlabel("Cycle");
  plt::ylabel("Error");
  plt::title("Goal-oriented error");
  plt::legend();
  plt::grid(true);

  // Save the plot
  plt::save(TEST_NAME+"-goal_oriented_convergence.png");
}

template <int dim>
void Problem<dim>::refine_mesh() {
  if(REFINEMENT_STRATEGY == "GO"){
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,error_indicators, 0.1, 0);
    
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