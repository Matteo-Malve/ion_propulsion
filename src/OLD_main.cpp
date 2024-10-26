#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h> // for the timer
#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out_faces.h>

#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>

#include <fstream>

using namespace dealii;
//using active_cell_iterator = typename DoFHandler<dim>::active_cell_iterator;
using std::cout;
using std::endl;


















template <int dim>
void Problem<dim>::SIMPLE_setup_system()
{
  primal_dof_handler.distribute_dofs(primal_fe);
 
  uh.reinit(primal_dof_handler.n_dofs());
  primal_rhs.reinit(primal_dof_handler.n_dofs());
 
  primal_constraints.clear(); 
  DoFTools::make_hanging_node_constraints(primal_dof_handler, primal_constraints);
  primal_constraints.close();
 
  DynamicSparsityPattern dsp(primal_dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(primal_dof_handler,
                                  dsp,
                                  primal_constraints,
                                  /*keep_constrained_dofs = */ false);
  primal_sparsity_pattern.copy_from(dsp);
  primal_system_matrix.reinit(primal_sparsity_pattern);
}

template <int dim>
void Problem<dim>::SIMPLE_assemble_system(){
  const QGauss<dim> quadrature_formula(/*primal_fe.degree + 1*/ 3);   // Aligned with primal solver
 
  FEValues<dim> fe_values(primal_fe,
                          quadrature_formula,
                          update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);
 
  const unsigned int dofs_per_cell = primal_fe.n_dofs_per_cell();
 
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
 
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
  for (const auto &cell : primal_dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
 
      cell_matrix = 0;
      cell_rhs    = 0;
 
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  (
                   fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                   fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                   fe_values.JxW(q_index));           // dx
 
              cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                              0.0 *                               // f(x)
                              fe_values.JxW(q_index));            // dx
            }
        }
 
      cell->get_dof_indices(local_dof_indices);
      primal_constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, primal_system_matrix, primal_rhs);
    }
		// Apply boundary values
  {
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values, ceiling_boundary_values;

    VectorTools::interpolate_boundary_values(primal_dof_handler,1, Functions::ConstantFunction<dim>(20000.), emitter_boundary_values);
    MatrixTools::apply_boundary_values(emitter_boundary_values, primal_system_matrix, uh, primal_rhs);

    VectorTools::interpolate_boundary_values(primal_dof_handler,2, Functions::ZeroFunction<dim>(), collector_boundary_values);
    MatrixTools::apply_boundary_values(collector_boundary_values, primal_system_matrix, uh, primal_rhs);

    // Ceiling
    // VectorTools::interpolate_boundary_values(primal_dof_handler,3, Functions::ZeroFunction<dim>(), ceiling_boundary_values);
    // MatrixTools::apply_boundary_values(ceiling_boundary_values, primal_system_matrix, uh, primal_rhs);
  }
}

template <int dim>
void Problem<dim>::SIMPLE_solve()
{
  SolverControl            solver_control(1000, 1e-6 * primal_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
 
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(primal_system_matrix, 1.2);
 
  solver.solve(primal_system_matrix, uh, primal_rhs, preconditioner);
 
  primal_constraints.distribute(uh);

  //cout<<"   Solved primal problem: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

template <int dim>
void Problem<dim>::SIMPLE_output_results() const
{
  ElectricFieldPostprocessor electric_field_postprocessor;
  DataOut<dim> data_out;
  data_out.attach_dof_handler(primal_dof_handler);
  data_out.add_data_vector(uh, "Potential");
  data_out.add_data_vector(uh, electric_field_postprocessor);
  data_out.build_patches();

	std::string filename;
  std::string meshName = extract_mesh_name();
	filename = std::string("automatic_lifting") + "-" + "primal" + "-" + meshName + "-" + Utilities::int_to_string(cycle, 1) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);
}

template <int dim>
void Problem<dim>::run()
{  
  // Create the mesh
  create_mesh(PATH_TO_MESH);

  // Cycles for refinement
	
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
      cout << "   Simple solver with Automatic imposition of BC " << " [FINAL]" << ':' << endl;
      SIMPLE_setup_system();
      SIMPLE_assemble_system();
      SIMPLE_solve();
      SIMPLE_output_results();
    } else {
      // Dual ---------------    
      cout<<"   Dual:"<<endl;
      setup_dual_system();
      assemble_dual_system();
      solve_dual();
      output_dual_results();
      cout<<"   Error estimation and Mesh refinement:"<<endl;
      refine_mesh();
    }
    ++cycle;
	}

  
}

int main(){
  try{
    Problem<2> iprop_problem;
    iprop_problem.run();
    }
  catch (std::exception &exc){
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

    return 1;
    }
  catch (...){
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
    }
 
  return 0;
}
