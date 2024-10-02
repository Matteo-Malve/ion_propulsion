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

// Physical Constants
const double eps_0 = 8.854; // 10^-12 [F/m]= [C^2 s^2 / kg / m^3]
const double eps_r = 1.0006;
const double Vmax = 2.e+4; // [V]
const double E_ON = 3.31e+6; // [V/m] corona inception threshold

const double R = 2.5e-4; // [m]
const double nn = 2;
const double L = nn*R; // [m]
const double X = 0.;//-L/2.; // [m]
const double g = 0.2; // [m]
const double mesh_height = 0.1;; // [m]

std::string PATH_TO_MESH = "../mesh/input_mesh.msh";
const unsigned int NUM_REFINEMENT_CYCLES = 3;
const unsigned int PRELIMINARY_REFINEMENT_STEPS = 0;

double get_emitter_height(const double &p)
{
	if (p <= X-L/2. || p >= X+L/2.)
		return 0.;

	double y = 0;
	double x = 0;

	const double left_center = X - L/2. + R;
	const double right_center = X + L/2. - R;

	if (p <= left_center)
		x = p - left_center;
	else if (p >= right_center)
		x = p - right_center;

	x /= R;
	y = R*std::sqrt( 1. - x * x);

	return y;
}

template <int dim>
class EmitterGeometry : public ChartManifold<dim, dim, dim-1>
{
public:
virtual Point<dim-1> pull_back(const Point<dim> &space_point) const override{		// space_point = p
	Point<dim-1> x;
	x[0] = space_point[0];
	if (dim == 3) 
		x[1] = space_point[2];
	return x;
}

virtual Point<dim> push_forward(const Point<dim-1> &chart_point) const override{		// chart_point = x
	const double y = get_emitter_height(chart_point[0]);
	Point<dim> p;
	p[0] = chart_point[0]; p[1] = y;
	if (dim == 3)
		p[2] = chart_point[1];
	return p;
}

virtual std::unique_ptr<Manifold<dim, dim>> clone() const override{
		return std::make_unique<EmitterGeometry<dim>>();
	}

};

template <int dim>
class RightHandSide : public Function<dim>
{
public:
	virtual void value_list( 	const std::vector< Point<dim>> &point_list, std::vector<double> &values, const unsigned int component = 0 ) const override {
	(void)component;
	AssertDimension (point_list.size(), values.size()); // Size check
		for (unsigned int p=0; p<point_list.size(); ++p)
			values[p] = 0.;
	}
};


template <int dim>
class Problem
{
public:
  Problem();

  void run();

private:
  void create_mesh(const std::string filename);

  void setup_primal_system();
  void setup_dual_system();
  void assemble_primal_system();
  void assemble_dual_system();
  void solve_primal();
  void solve_dual();
  void output_primal_results(const unsigned int cycle);
  void output_dual_results(const unsigned int cycle);

  double estimate_error(Vector<float> &error_indicators) const;
  void refine_mesh();

	void SIMPLE_setup_system();
  void SIMPLE_assemble_system();
  void SIMPLE_solve();
  void SIMPLE_output_results(const unsigned int cycle) const;

  Triangulation<dim> triangulation;

  FE_Q<dim>       primal_fe;
  FE_Q<dim>       dual_fe;

  QGauss<dim> dual_quadrature;

	DoFHandler<dim> primal_dof_handler;
  DoFHandler<dim> dual_dof_handler;

  AffineConstraints<double> primal_constraints; // to deal with hanging nodes
  SparsityPattern      primal_sparsity_pattern;
  SparseMatrix<double> primal_system_matrix;

  Vector<double> primal_rhs;
  Vector<double> uh0;
	Vector<double> Rg_vector;
  Vector<double> primal_solution;


  AffineConstraints<double> dual_constraints; // to deal with hanging nodes

  SparsityPattern      dual_sparsity_pattern;
  SparseMatrix<double> dual_system_matrix;

  Vector<double> dual_solution;
  Vector<double> dual_rhs;

  RightHandSide<dim> rhs_function; // At the moment, this is equivalent to Functions::ZeroFunction<dim>()

  Timer timer;

  class ElectricFieldPostprocessor;
  class IonizationAreaPostprocessor;
};

template <int dim>
class Problem<dim>::ElectricFieldPostprocessor : public DataPostprocessorVector<dim>
{
public:
  ElectricFieldPostprocessor ():  DataPostprocessorVector<dim> ("Electric_Field", update_gradients) {}

  virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
									 std::vector<Vector<double> > &computed_quantities) const override 
    {
    AssertDimension (input_data.solution_gradients.size(), computed_quantities.size()); // size check
    for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p){
      AssertDimension (computed_quantities[p].size(), dim); // dimension check
      for (unsigned int d=0; d<dim; ++d)
        computed_quantities[p][d] = -input_data.solution_gradients[p][d];
      }
    }
};

template <int dim>
class Problem<dim>::IonizationAreaPostprocessor : public DataPostprocessorScalar<dim>
{
public:
  IonizationAreaPostprocessor ():  DataPostprocessorScalar<dim> ("normalized_overfield", update_gradients) {}

  virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
									 std::vector<Vector<double>> &computed_quantities) const override {
	  AssertDimension (input_data.solution_values.size(), computed_quantities.size()); // Size check

	  for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
	  {
	AssertDimension(computed_quantities[p].size(), 1);
		  double magnitude = 0;
		  for (unsigned int d=0; d<dim; ++d)
			  magnitude += input_data.solution_gradients[p][d] * input_data.solution_gradients[p][d];

		  const double overfield = std::max(0., std::sqrt(magnitude) - E_ON);

		  computed_quantities[p](0) = overfield/E_ON;
	  }
  }
};


template <int dim>
class Evaluate_Rg : public Function<dim>
{
public:
  virtual double value(const Point<dim> & p,
					   const unsigned int component = 0) const override{
	  (void)component;

	  const auto y = p[1];
	  const auto x = p[0];

	if (x <= (X-L/2.-2.*R) || x >= (X+L/2.+2.*R) || y >= 3.*R)
		return 0.;

	double Rg = 0.;
	double d = 0.;
	const double left_center = X - L/2. + R;
	const double right_center = X + L/2. - R;

	if (x >= right_center)
		d = sqrt((x-right_center)*(x-right_center) + y*y);
	else if (x < right_center && x > left_center )
		d = y;
	else if (x <= left_center )
		d = sqrt((x-left_center)*(x-left_center) + y*y);
	else
		cout << "ERROR! Not implemented!" << endl;

	d = std::max(d,R); /* As the boundary only approximates a circle,
	* some points on the emitter edge will be inside the radius,
	* and would produce an Rg higher than Vmax with the current function choice */

	if (d < 3.*R)
		Rg = Vmax * (1. - (d-R)/(2.*R)) * (1. - (d-R)/(2.*R)); //* (1. - (d-R)/(2.*R))

		  if (Rg > Vmax) cout << "Error! Rg is " << Rg << " in " << x << ", " << y << endl;

		  return Rg;
	  }
};

Tensor<1,2> emitter_normal(const Point<2> p) {
	//Compute the normal at point p of a circular emitter of length L
	// with circular edges of radius R centered in [-R,0] and in [-L+R,0]

	Tensor<1,2> normal;
	normal[0] = 0.;
	normal[1] = p[1];

	const double x = p[0];
	const double left_center = X - L/2. + R;
	const double right_center = X + L/2. - R;

	if (x >= right_center)
		normal[0] = x - right_center;
	else if (x <= left_center)
		normal[0] = x - left_center;

	const double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]); // The norm of this vector should always be nearly equal to R

	return normal/norm;
}


static auto evaluate_grad_Rg = [](const double x, const double y) {

	Tensor<1,2> grad_Rg;

	// Gradient computed analytically by hand. Check Evaluate_Rg.h for the primitive function.
	grad_Rg[0] = 0.;
	grad_Rg[1] = 0.;

	if (x <= (X-L/2.-2.*R) || x >= (X+L/2.+2.*R) || y >= 3.*R)
		return grad_Rg;

	double d = 0.;
	const double left_center = X - L/2. + R;
	const double right_center = X + L/2. - R;

	if (x >= right_center )
		d = std::sqrt((x-right_center)*(x-right_center) + y*y);
	else if (x < right_center && x > left_center )
		d = y;
	else if (x <= left_center )
		d = std::sqrt((x-left_center)*(x-left_center) + y*y);
	else
		cout << "ERROR! Not implemented!" << endl;

	d = std::max(d,R);

	if (d < 3.*R) {
		if (x >= right_center ) {
			grad_Rg[0] = - Vmax/R * (1. - (d-R)/(2.*R)) * (x-right_center)/d; //* (1. - (d-R)/(2.*R)) * 1.5
			grad_Rg[1] = - Vmax/R * (1. - (d-R)/(2.*R)) * y/d;
		} else if (x < right_center && x > left_center ) {
			grad_Rg[0] = 0.;
			grad_Rg[1] = - Vmax/R * (1. - (d-R)/(2.*R));
		} else if (x <= left_center ) {
			grad_Rg[0] = - Vmax/R * (1. - (d-R)/(2.*R)) * (x-left_center)/d;
			grad_Rg[1] = - Vmax/R * (1. - (d-R)/(2.*R)) * y/d;
		}
	}

	  return grad_Rg;
};


template <int dim>
Problem<dim>::Problem()
  : primal_fe(1)
  , dual_fe(2)
  , dual_quadrature(dual_fe.degree + 1)
  , primal_dof_handler(triangulation)
  , dual_dof_handler(triangulation)
  //, mapping()
{}

template <int dim>
void Problem<dim>::create_mesh(const std::string filename)
{
	cout << endl << "Reading file: " << filename << endl;
	std::ifstream input_file(filename);
	GridIn<2>       grid_in;
	grid_in.attach_triangulation(triangulation);
	grid_in.read_msh(input_file);
  
	const types::manifold_id emitter = 1;
	EmitterGeometry<2> emitter_manifold;
  
	for (auto &cell : triangulation.active_cell_iterators()) {
		if (cell->at_boundary()) {
			for (unsigned int f=0; f<4; ++f) {
				if (cell->face(f)->at_boundary()) {

					const Point<dim> fc = cell->face(f)->center();

					if (fc[1] > 0 && fc[1] < 3*R && std::abs(fc[0]-X) < L/2) {
						for (unsigned int v = 0; v < 2; ++v)
							cell->face(f)->vertex(v)[1] =  get_emitter_height(cell->face(f)->vertex(v)[0]);
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
	
	for (unsigned int i = 0; i < PRELIMINARY_REFINEMENT_STEPS; ++i) {

		Vector<float> criteria(triangulation.n_active_cells());
		cout  << "Active cells " << triangulation.n_active_cells() << endl;
		unsigned int ctr = 0;

		for (auto &cell : triangulation.active_cell_iterators()) {

			const Point<dim> c = cell->center();
			const double d = std::sqrt( (c[0]-X)*(c[0]-X) + c[1]*c[1]);

			if ( d <= PRELIMINARY_REFINEMENT_STEPS/(i+1)*2.*R)
				criteria[ctr++] = 1;
			else
				criteria[ctr++] = 0;
		}
		GridRefinement::refine(triangulation, criteria, 0.5);
		triangulation.execute_coarsening_and_refinement();
		cout<<"   Executed preliminary coarsening and refinement"<<endl;
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

  /*
	cout  << "Final number of active cells: " << triangulation.n_active_cells() << endl;
	cout <<"   Executed Pre-Refinement"<< endl;
  */
}

template <int dim>
void Problem<dim>::setup_primal_system()
{
  primal_dof_handler.distribute_dofs(primal_fe);
  
	Rg_vector.reinit(primal_dof_handler.n_dofs());
	uh0.reinit(primal_dof_handler.n_dofs());
  primal_solution.reinit(primal_dof_handler.n_dofs());
  primal_rhs.reinit(primal_dof_handler.n_dofs());
 
  primal_constraints.clear();
	DoFTools::make_hanging_node_constraints(primal_dof_handler, primal_constraints);
	primal_constraints.close();


  DynamicSparsityPattern dsp(primal_dof_handler.n_dofs(),primal_dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(primal_dof_handler, dsp);
  primal_constraints.condense(dsp);
  /*
  DynamicSparsityPattern dsp(primal_dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(primal_dof_handler, dsp, primal_constraints, false);
  */
  primal_sparsity_pattern.copy_from(dsp);
  primal_system_matrix.reinit(primal_sparsity_pattern);

  std::cout << "      Number of degrees of freedom: " << primal_dof_handler.n_dofs()<<std::endl;
  
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

  std::cout << "      Number of degrees of freedom: " << dual_dof_handler.n_dofs()<< std::endl;
}


template <int dim>
void Problem<dim>::assemble_primal_system()
{
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
    auto grad_Rg_xq = evaluate_grad_Rg(quad_point_coords[0],quad_point_coords[1]);

    // assemble A_loc(Rg,v)
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      cell_rhs(i) -= eps_r*eps_0*
              (Rg_fe_values.shape_grad(i, q_point) *     // grad phi_i(x_q)
              grad_Rg_xq *                               // grad_Rg(x_q)
              Rg_fe_values.JxW(q_point));                // dx
  }
      // Local to global
      cell->get_dof_indices(local_dof_indices);
      //primal_constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, primal_system_matrix, primal_rhs);
      for (const unsigned int i : fe_values.dof_indices()){
        for (const unsigned int j : fe_values.dof_indices())
          primal_system_matrix.add(local_dof_indices[i], 
                            local_dof_indices[j],
                            cell_matrix(i, j)); 
        primal_rhs(local_dof_indices[i]) += cell_rhs(i);
      }
  }

  // Apply boundary values
  {
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    VectorTools::interpolate_boundary_values(primal_dof_handler,1, Functions::ZeroFunction<dim>(), emitter_boundary_values);
    MatrixTools::apply_boundary_values(emitter_boundary_values, primal_system_matrix, primal_solution, primal_rhs);

    VectorTools::interpolate_boundary_values(primal_dof_handler,2, Functions::ZeroFunction<dim>(), collector_boundary_values);
    MatrixTools::apply_boundary_values(collector_boundary_values, primal_system_matrix, primal_solution, primal_rhs);
  }
  primal_constraints.condense(primal_system_matrix);
  primal_constraints.condense(primal_rhs);

}

template <int dim>
void Problem<dim>::assemble_dual_system()
{
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
        Tensor<1, dim> n = emitter_normal(face->center());

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
}

template <int dim>
void Problem<dim>::solve_primal()
{
  // Solver setup
  const unsigned int it_max = 1e+4;
  const double rel_tol = 1.e-6*primal_rhs.l2_norm();
  const double abs_tol = 1.e-12;
  const double tol = abs_tol + rel_tol;
  SolverControl            solver_control(it_max, tol);
  SolverCG<Vector<double>> solver(solver_control);

  // Solve linear system
  double relaxation_parameter = 1.2;
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(primal_system_matrix, relaxation_parameter);
  solver.solve(primal_system_matrix, primal_solution, primal_rhs, preconditioner);

  // uh0
	uh0 = primal_solution;
	primal_constraints.distribute(uh0);
  primal_constraints.distribute(primal_solution);

	// Retrieve lifting
  VectorTools::interpolate(primal_dof_handler,Evaluate_Rg<dim>(),Rg_vector);
  primal_solution += Rg_vector;               // uh = u0 + Rg

  cout<<"      Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
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
  cout<<"      Solved system: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

std::string extract_mesh_name() {
  // Find the last '/' character to isolate the filename
  size_t lastSlash = PATH_TO_MESH.rfind('/');
  // Find the last '.' character to remove the file extension
  size_t lastDot = PATH_TO_MESH.rfind('.');
  // Extract the substring between the last '/' and the last '.'
  std::string meshName = PATH_TO_MESH.substr(lastSlash + 1, lastDot - lastSlash - 1);
  return meshName;
}

template <int dim>
void Problem<dim>::output_primal_results(const unsigned int cycle)
{
  // const Point<dim> evaluation_point(0.5*g, 0.1*g);
  // const double x = VectorTools::point_value(primal_dof_handler, primal_solution, evaluation_point);
  // cout << "   Potential at sample point (" << evaluation_point[0] << "," << evaluation_point[1] << "): " << x << endl;

	ElectricFieldPostprocessor electric_field_postprocessor;
	IonizationAreaPostprocessor ionization_area_postprocessor;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(primal_dof_handler);
  data_out.add_data_vector(primal_solution, "Potential");
	data_out.add_data_vector(uh0, "uh0");
	data_out.add_data_vector(Rg_vector, "Rg");
  data_out.add_data_vector(primal_solution, electric_field_postprocessor);
	data_out.add_data_vector(primal_solution, ionization_area_postprocessor);
  data_out.build_patches(); // mapping

  std::string filename;
  std::string meshName = extract_mesh_name();
	filename = std::string("Rg_manual") + "-" + "primal" + "-" + meshName + "-" + Utilities::int_to_string(cycle, 1) + ".vtk";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtk(output);

  /* Output the mesh in .msh format (Gmsh)
  std::ofstream out(meshName+"-"+Utilities::int_to_string(cycle, 1)+".msh");
  GridOut grid_out;
  grid_out.write_msh(triangulation, out);

  std::cout << "Mesh saved to mesh.msh" << std::endl;
  */
}

template <int dim>
void Problem<dim>::output_dual_results(const unsigned int cycle)
{
	ElectricFieldPostprocessor electric_field_postprocessor;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dual_dof_handler);
  data_out.add_data_vector(dual_solution, "Potential");
  data_out.add_data_vector(dual_solution, electric_field_postprocessor);
  data_out.build_patches(); // mapping

  std::string filename;
  std::string meshName = extract_mesh_name();

	filename =  std::string("Rg_manual") + "-" + "dual-" + meshName + "-" + Utilities::int_to_string(cycle, 1) + ".vtk";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtk(output);
}


template <int dim>
void Problem<dim>::refine_mesh(){
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

	GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,error_indicators, 0.8, 0.0); // CHANGED FROM: 0.8, 0.02

  // Execute refinement
  triangulation.execute_coarsening_and_refinement();
  // cout<<"      Executed coarsening and refinement"<<endl;  // redundant message
}

template <int dim>
double Problem<dim>::estimate_error(Vector<float> &error_indicators) const
{
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
    VectorTools::interpolate(dual_dof_handler,Evaluate_Rg<dim>(),Rg_vector);

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
            auto grad_Rg_xq = evaluate_grad_Rg(quad_point_coords[0],quad_point_coords[1]);
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
void Problem<dim>::SIMPLE_setup_system()
{
  primal_dof_handler.distribute_dofs(primal_fe);
 
  primal_solution.reinit(primal_dof_handler.n_dofs());
  primal_rhs.reinit(primal_dof_handler.n_dofs());
 
  primal_constraints.clear();   // clear from previous cycle

  DoFTools::make_hanging_node_constraints(primal_dof_handler, primal_constraints);
	/* Alternative: Impose BCs as AffineConstraints rather than mapping
	
  VectorTools::interpolate_boundary_values(primal_dof_handler,
                                           types::boundary_id(1),
                                           Functions::ConstantFunction<dim>(20000.),
                                           primal_constraints);
  VectorTools::interpolate_boundary_values(primal_dof_handler,
                                           types::boundary_id(2),
                                           Functions::ZeroFunction<dim>(),
                                           primal_constraints);               */                          
 
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
  const QGauss<dim> quadrature_formula(primal_fe.degree + 1);
 
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
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    VectorTools::interpolate_boundary_values(primal_dof_handler,1, Functions::ConstantFunction<dim>(20000.), emitter_boundary_values);
    MatrixTools::apply_boundary_values(emitter_boundary_values, primal_system_matrix, primal_solution, primal_rhs);

    VectorTools::interpolate_boundary_values(primal_dof_handler,2, Functions::ZeroFunction<dim>(), collector_boundary_values);
    MatrixTools::apply_boundary_values(collector_boundary_values, primal_system_matrix, primal_solution, primal_rhs);
  }
}

template <int dim>
void Problem<dim>::SIMPLE_solve()
{
  SolverControl            solver_control(1000, 1e-6 * primal_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
 
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(primal_system_matrix, 1.2);
 
  solver.solve(primal_system_matrix, primal_solution, primal_rhs, preconditioner);
 
  primal_constraints.distribute(primal_solution);

  //cout<<"   Solved primal problem: "<<solver_control.last_step()  <<" CG iterations needed to obtain convergence." <<endl;
}

template <int dim>
void Problem<dim>::SIMPLE_output_results(const unsigned int cycle) const
{
  ElectricFieldPostprocessor electric_field_postprocessor;
  DataOut<dim> data_out;
  data_out.attach_dof_handler(primal_dof_handler);
  data_out.add_data_vector(primal_solution, "Potential");
  data_out.add_data_vector(primal_solution, electric_field_postprocessor);
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
	unsigned int cycle = 0;
	while (cycle < NUM_REFINEMENT_CYCLES) {
    cout << endl << "Cycle " << cycle << ':' << endl;
    std::cout << "   Number of active cells:       "<< triangulation.n_active_cells() << std::endl;

    // Primal --------------
		cout<<"   Primal:"<<endl;
		setup_primal_system();
		assemble_primal_system();
		solve_primal();
		output_primal_results(cycle);

    // Dual ---------------    
		cout<<"   Dual:"<<endl;
		setup_dual_system();
		assemble_dual_system();
		solve_dual();
		output_dual_results(cycle);

    // --------------------  
		cout<<"   Error estimation and Mesh refinement:"<<endl;
		refine_mesh();

		++cycle;
	}

  // Cycles for solution
  cout << "Cycle " << cycle << " [FINAL]" << ':' << endl;
	SIMPLE_setup_system();
	SIMPLE_assemble_system();
	SIMPLE_solve();
	SIMPLE_output_results(cycle);
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
