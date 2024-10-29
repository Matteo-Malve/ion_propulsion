#ifndef PROBLEM_H
#define PROBLEM_H

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

#include <fstream>




#include "dual_functional.h"
#include "globals.h"
#include "postprocessors.h"
#include "utilities.h"

using namespace dealii;

template <int dim>
class Problem {
public:
	Problem();
	void run();

private:
	void create_mesh();
	void setup_primal_system();
	void assemble_primal_system();
	void solve_primal();
	void output_primal_results();

	void setup_dual_system();
	void assemble_dual_system();
	void solve_dual();
	void output_dual_results();

	void estimate_error();
	void refine_mesh();

	void test_convergence();

	Triangulation<dim>            triangulation;

	DoFHandler<dim>               primal_dof_handler, 
																dual_dof_handler;
																
	FE_Q<dim>                     primal_fe, 
																dual_fe;

	AffineConstraints<double> 		primal_constraints, 
																dual_constraints;
	
	SparsityPattern								primal_sparsity_pattern,
																dual_sparsity_pattern;
	
	SparseMatrix<double> 					primal_system_matrix, 
																dual_system_matrix;

	Vector<double> 								primal_rhs, 
																dual_rhs;

	Vector<double> 								uh0, 
																uh, 
																zh;

	Vector<double> 								Rg_primal, 
																Rg_dual;

	Vector<double> 								uh0_on_dual_space,
																Rg_plus_uh0hat, 
																error_indicators;															
	
	RightHandSide4<dim> 					rhs_function;
	ExactSolution4<dim>						exact_solution_function;

	int cycle;
};

#endif
