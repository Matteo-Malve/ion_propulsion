/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2021 by the deal.II authors
 *
 * This file is based on step-6 of the examples section of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Matteo Menessini, Politecnico di Milano, 2023
 *
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h> // for the timer
#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h> // To use manifolds
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

#include <fstream>


using namespace dealii;

// Dichiarazione variabili e funzioni

template <int dim>
class Problem
{
public:
	enum EmitterType
	{
		wire,
		blade,
		plate,
		tip
	};

  Problem(const EmitterType emitter_type);

  void run();

private:
  void create_mesh();
  void setup_system();
  void solve();
  void refine_grid(const unsigned int min_grid_level,
          	  	  const unsigned int max_grid_level);
  void output_results();

  Triangulation<dim> triangulation;

  FE_Q<dim>       fe;
  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints; // to deal with hanging nodes

  SparseMatrix<double> system_matrix;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> laplace_matrix;
  //SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
  //Vector<double> old_solution;

  Timer timer;

  const EmitterType emitter_type;

  unsigned int cycle = 0;
  const unsigned int Nmax = 20; // maximum number of refinement cycles

  Vector<float> values;
  const float conv_tol = 1e-3;
};

// Qui alcune variabili globali usate in diversi punti del codice
// (non è il massimo dichiararle così...)

// Wire
const double wire_radius = 0.0025; // [cm]

// Plate & Blade
const double emitter_height = 0.02;// [cm]
const double emitter_length = 0.2;// [cm]

// Blade
const double blade_angle = 25. / 180 * numbers::PI; // [rad] typically 22°-30°
const double t_blade = tan(blade_angle);
const double blade_length = emitter_height/t_blade;

const double eps0 = 8.854*1e-12; // [F/m]
const double electrode_distance = 2.;
const double airfoil_length = 10.;
const double collector_length = 3.;


// Funzione per determinare l'altezza del profilo NACA al variare dell'ascissa
double get_height(const double &x)
{
	double y = 0;

	if ( std::fabs(x-1) > 1e-12 && std::fabs(x) > 1e-12 ) {
		// 4-Digit NACA:
		double a0 = 0.2969;
		double a1 = -0.126;
		double a2 = -0.3516;
		double a3 = 0.2843;
		double a4 = -0.1036;
		double t = 0.5; // Last 2 digits of the NACA divided by 20

		y = t*( a0 * std::sqrt(x) + a1 * x + a2 * pow(x,2.0) + a3 * pow(x,3.0) + a4 * pow(x,4.0) );
	}

	return y;
}


// Classe che definisce il "manifold" del NACA
template <int dim, int sign>
class AirfoilGeometry : public ChartManifold<dim, dim, dim-1>
{
public:
	virtual Point<dim-1> pull_back(const Point<dim> &space_point) const override;

	virtual Point<dim> push_forward(const Point<dim-1> &chart_point) const override;

	virtual std::unique_ptr<Manifold<dim, dim>> clone() const override;

};

template <int dim, int sign>
std::unique_ptr<Manifold<dim, dim>> AirfoilGeometry<dim, sign>::clone() const
{
  return std::make_unique<AirfoilGeometry<dim, sign>>();
}

template <int dim, int sign>
Point<dim> AirfoilGeometry<dim, sign>::push_forward(const Point<dim-1>  &x) const
{
	const double y = sign*get_height( (x[0] - electrode_distance)/airfoil_length )*airfoil_length;

	Point<dim> p;
	p[0] = x[0]; p[1] = y;

	if (dim == 3) {
		p[2] = x[1];
	}

	return p;
}

template <int dim, int sign>
Point<dim-1>  AirfoilGeometry<dim, sign>::pull_back(const Point<dim> &p) const
{
	Point<dim-1> x;
	x[0] = p[0];

	if (dim == 3) {
		x[1] = p[2];
	}

	return x;
}


//Classe che definisce il "manifold" della lama
template <int dim, int sign>
class BladeGeometry : public ChartManifold<dim, dim, dim-1>
{
public:
	virtual Point<dim-1> pull_back(const Point<dim> &space_point) const override;

	virtual Point<dim> push_forward(const Point<dim-1> &chart_point) const override;

	virtual std::unique_ptr<Manifold<dim, dim>> clone() const override;

};

template <int dim, int sign>
std::unique_ptr<Manifold<dim, dim>> BladeGeometry<dim, sign>::clone() const
{
  return std::make_unique<BladeGeometry<dim, sign>>();
}

template <int dim, int sign>
Point<dim> BladeGeometry<dim, sign>::push_forward(const Point<dim-1>  &x) const
{
	double y = 0;

	if (sign == 2)
		y = std::max(-emitter_height, t_blade*x[0]);
	else
		y = sign*std::min(emitter_height/2, -t_blade*x[0]);

	Point<dim> p;
	p[0] = x[0]; p[1] = y;

	if (dim == 3) {
		p[2] = x[1];
	}

	return p;
}

template <int dim, int sign>
Point<dim-1>  BladeGeometry<dim, sign>::pull_back(const Point<dim> &p) const
{
	Point<dim-1> x;
	x[0] = p[0];

	if (dim == 3) {
		x[1] = p[2];
	}

	return x;
}

//  Funzione per creare la mesh del profilo NACA
void CreateNACA(const double mesh_height, const double airfoil_length,
				 const double electrode_distance, Triangulation<2> &mesh)
{
// Input data
	const double cube_side = electrode_distance/2;
	const double left_pad = mesh_height/2 - electrode_distance/2;
	const double tail_length = left_pad - airfoil_length  - electrode_distance/2;

	Triangulation<2> tria, naca_up, naca_down, aux;

	// NACA triangulation
	const Point<2> top_left(cube_side, mesh_height/2);
	const Point<2> middle(electrode_distance+airfoil_length+tail_length, 0.);
	const Point<2> bottom_left(cube_side, -mesh_height/2);

	unsigned int a = (int) middle(0) - top_left(0);
	unsigned int b = (int) top_left(1) - middle(1);

	GridGenerator::subdivided_hyper_rectangle( naca_up, {a,b}, top_left, middle);
	GridGenerator::subdivided_hyper_rectangle( naca_down, {a,b}, middle, bottom_left);

	const double left_limit = electrode_distance;
	const double right_limit = electrode_distance + airfoil_length;
	const double x_unit = (middle(0) - top_left(0))/((double) a);

	  for (auto &face : naca_up.active_face_iterators())
	  {
		  const Point<2> c = face->center();
		  if ( (c(1) == (int) c(1)) ) {
		  if (face->at_boundary())
		  {
			  if ( (c(0) > left_limit) && (c(0) < right_limit) && (c(1) < airfoil_length/10)) //
			  {
					for (const auto i : face->vertex_indices())
					{
					  Point<2> &v = face->vertex(i);
					  double h = get_height( (v(0) - electrode_distance)/airfoil_length );
					  v(1) = h * airfoil_length;
					}
			  }
		  } else { // if !face->at_boundary()
			  if ( (c(0) > left_limit - x_unit) && (c(0) < right_limit + x_unit ))
			  {
				for (const auto i : face->vertex_indices())
				{
				  Point<2> &v = face->vertex(i);
				  double h = get_height( (v(0) - electrode_distance + x_unit)/(airfoil_length+2*x_unit) );
				  v(1) += h * airfoil_length * (mesh_height-c(1)*2)/mesh_height;
				}
			  }
		  }
	  }
	 }

    for (auto &face : naca_down.active_face_iterators()) {
		  const Point<2> c = face->center();
		  if ( (c(1) == (int) c(1)) ) {
			  if (face->at_boundary())
			  {
				  if ( (c(0) > left_limit) && (c(0) < right_limit) && (c(1) > - airfoil_length/10))
				  {
						for (const auto i : face->vertex_indices())
						{
						  Point<2> &v = face->vertex(i);
						  double h = get_height( (v(0) - electrode_distance)/airfoil_length );
						  v(1) = -h * airfoil_length;
						}
				  }
			  } else { // if !face->at_boundary()
				  if ( (c(0) > left_limit - x_unit) && (c(0) < right_limit + x_unit ))
				  {
					for (const auto i : face->vertex_indices())
					{
					  Point<2> &v = face->vertex(i);
					  double h = get_height( (v(0) - electrode_distance + x_unit)/(airfoil_length+2*x_unit) );
					  v(1) -= h * airfoil_length * (mesh_height+c(1)*2)/mesh_height;
					}
				  }
			  }
		  }
	  }

	GridGenerator::merge_triangulations( naca_down, naca_up, mesh);
}


// Funzione per assegnare i vari manifold e il numero identificativo dei boundaries
template <int dim>
void SetManifoldsAndBoundaries(Triangulation<dim> &mesh)
{
    Assert(dim == 2, ExcNotImplemented());

	const types::manifold_id  	airfoil_id_top = 11;
	const types::manifold_id  	airfoil_id_bottom = 12;

	AirfoilGeometry<dim,1> airfoil_manifold_top;
	mesh.set_manifold(airfoil_id_top, airfoil_manifold_top);

	AirfoilGeometry<dim,-1> airfoil_manifold_bottom;
	mesh.set_manifold(airfoil_id_bottom, airfoil_manifold_bottom);

	//Set boundaries
	const types::boundary_id airfoil_id = 10;
	const types::boundary_id collector_id = 2;
	const types::boundary_id emitter_id = 1;

	  for (auto &face : mesh.active_face_iterators())
	  {
		  if (face->at_boundary())
		  {
			  const Point<dim> c = face->center();

			  if ( (c(0) > electrode_distance) && (c(0) < electrode_distance + airfoil_length)
					  && (std::fabs(c(1)) <= airfoil_length/10 + 1e-12))
			  {
				  face->set_boundary_id(airfoil_id);

				  if (c(0) < electrode_distance + collector_length)
					  face->set_boundary_id(collector_id);

				  if (c(1) > 0)
						face->set_manifold_id(airfoil_id_top);
				  else
						face->set_manifold_id(airfoil_id_bottom);
			  }

			  if ( c.square() <= electrode_distance/2)
				  face->set_boundary_id(emitter_id);
		  }
	  }

	  const Point<2>             center(electrode_distance + airfoil_length/10, 0.);
	  const SphericalManifold<2> shell(center);

	  const types::manifold_id     shell_id = 8;
	  mesh.set_manifold(shell_id, shell);

	  for (auto &cell : mesh.active_cell_iterators())
	  {
		  if (cell->at_boundary())
		  {
			  const Point<dim> c = cell->center();

			  if ( (c(0) > electrode_distance) && (c(0) < electrode_distance + airfoil_length/10) && (std::fabs(c(1)) < airfoil_length/10))
			  {

				  for (unsigned int i = 0; i < 4; ++i) {
					  if ( !(cell->face(i)->at_boundary()) )
						  cell->face(i)->set_manifold_id(shell_id);
				  }
			  }
		  }
	  }
}

// Funzione per calcolo norma di tensore
template <int dim>
double L2Norm(const Tensor<1,dim> &input)
{
	double magnitude = 0;
	for (unsigned int i=0; i<dim; ++i)
		magnitude += input[i]*input[i];

	return std::sqrt(magnitude);
}


// Funzione usata in post-processing per estrarre il campo elettrico come -gradiente del potenziale
template <int dim>
class GradientPostprocessor : public DataPostprocessorVector<dim>
{
public:
  GradientPostprocessor ():  DataPostprocessorVector<dim> ("E", update_gradients) {}

  virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
		  	  	  	  	  	  	  	 std::vector<Vector<double> > &computed_quantities) const override {
    AssertDimension (input_data.solution_gradients.size(), computed_quantities.size()); // size check

    for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
      {
        AssertDimension (computed_quantities[p].size(), dim); // dimension check
        for (unsigned int d=0; d<dim; ++d)
          computed_quantities[p][d] = -input_data.solution_gradients[p][d];
      }
  }
};

// !!! Scelta degli elementi finiti: lineari o quadratici
template <int dim>
Problem<dim>::Problem(const EmitterType emitter_type)
  : fe(1) // linear (1) or quadratic (2) elements
  , dof_handler(triangulation)
  , emitter_type(emitter_type)
{}


// Creazione della griglia computazionale
template <int dim>
void Problem<dim>::create_mesh()
{
	const double mesh_height = 50.; // Deve essere pari
	const double naca_length = 10.;

	Triangulation<2> naca, emitter_mesh;
	CreateNACA(mesh_height, naca_length, electrode_distance, naca);
	const double cube_side = electrode_distance/2;
	const double left_pad = mesh_height/2 - electrode_distance/2;

    std::cout << std::endl;


    // Switch per i vari tipi di emettitori
    switch (emitter_type)
      {
        case wire:
          {
        	std::cout << "Running for WIRE emitter with radius: " << wire_radius*10 << " mm" << std::endl;
        	GridGenerator::plate_with_a_hole( emitter_mesh, wire_radius, cube_side, mesh_height/2-cube_side,
        									mesh_height/2-cube_side, left_pad, 0, {0.,0.});
            break;
          }

        case blade:
                  {
                  	const double th = emitter_height*1e-6;
          			const double radius = cube_side/3;

                  	std::cout << "Running for BLADE emitter with length: " << emitter_length*10  << "mm, height: "
                  			<< emitter_height*10 << "mm, angle: " << blade_angle*180/numbers::PI << "°" << std::endl;
          			GridGenerator::plate_with_a_hole( emitter_mesh, radius, cube_side, mesh_height/2-cube_side,
          											mesh_height/2-cube_side, left_pad, 0, {0.,0.});
          			emitter_mesh.reset_all_manifolds();
          			emitter_mesh.set_all_manifold_ids_on_boundary(0);

          			for (auto &face : emitter_mesh.active_face_iterators())
          			{
          				if (face->at_boundary()) {
          				  const Point<2> c = face->center();

          				  if (std::sqrt(c.square()) <= cube_side) {
          					for (const auto i : face->vertex_indices())
          					{
          					  Point<2> &v = face->vertex(i);
          					  if ( std::fabs(std::sqrt(v.square()) - radius) < th) {
          						  if (std::fabs(v(1)) < th) {
          							  v(1) = 0.-emitter_height/2*(v(0) < th);
          							  v(0) = 0. - emitter_length*(v(0) < th);
          						  } else {
          							  if (std::fabs(v(0)) < th) {
          								  v(0) = -2*emitter_length/3*(v(1) > th) -(emitter_length+blade_length)/2*(v(1) < th);
          								  v(1) = 0. - emitter_height*(v(1) < th);
          							  } else {
          								  v(0) =  -blade_length*(v(0) > th)*(v(1) < th) -emitter_length*(v(0) < th) -emitter_length/3*(v(0) > th)*(v(1) > th);
          								  v(1) = 0. - emitter_height*(v(1) < th);
          							  }

          						  }
          					   }
          					 }
          				  }
          				}
          			}
                    break;
                  }

        case plate:
		  {
		    const double th = emitter_height*1e-6;

        	std::cout << "Running for PLATE emitter with length: " << emitter_length*10  << "mm, height: "
        			<< emitter_height*10 << "mm" << std::endl;

			const double radius = cube_side/3;
			GridGenerator::plate_with_a_hole( emitter_mesh, radius, cube_side, mesh_height/2-cube_side,
											mesh_height/2-cube_side, left_pad, 0, {0.,0.});
			emitter_mesh.reset_all_manifolds();
			emitter_mesh.set_all_manifold_ids_on_boundary(0);

			for (auto &face : emitter_mesh.active_face_iterators())
				{
					if (face->at_boundary()) {
					  const Point<2> c = face->center();

					  if (std::sqrt(c.square()) <= cube_side) {
						for (const auto i : face->vertex_indices())
						{
						  Point<2> &v = face->vertex(i);

						  if (std::fabs(std::sqrt(v.square()) - radius) < th)
						  {
							  if (std::fabs(v(1)) > th) {
								  v(1) = emitter_height/2 * ((v(1) > 0) - (v(1) < 0));
								  v(0) = 0 * (v(0) > th) - emitter_length * (v(0) < -th) - emitter_length/2 * (std::fabs(v(0)) <= th);
							  } else
								  v(0) = 0. - emitter_length*(v(0) < 0);
						  }

						}

					  }
					}
				}
			break;
		  }

        case tip:
		  {
		    const double th = emitter_height*1e-6;

        	std::cout << "Running for TIP emitter with length: " << emitter_length*10  << "mm, height: "
        			<< emitter_height*10 << "mm, angle: " << blade_angle*360/numbers::PI << "°" << std::endl;

			const double radius = cube_side/3;
			GridGenerator::plate_with_a_hole( emitter_mesh, radius, cube_side, mesh_height/2-cube_side,
											mesh_height/2-cube_side, left_pad, 0, {0.,0.});
			emitter_mesh.reset_all_manifolds();
			emitter_mesh.set_all_manifold_ids_on_boundary(0);

			for (auto &face : emitter_mesh.active_face_iterators())
				{
					if (face->at_boundary()) {
					  const Point<2> c = face->center();

					  if (std::sqrt(c.square()) <= cube_side) {
						for (const auto i : face->vertex_indices())
						{
						  Point<2> &v = face->vertex(i);

						  if (std::fabs(std::sqrt(v.square()) - radius) < th)
						  {
							  if (std::fabs(v(1)) > th) {
								  v(1) = emitter_height/2 * ((v(1) > 0) - (v(1) < 0));
								  v(0) = -blade_length * (v(0) > th) - emitter_length * (v(0) < -th) - emitter_length/2 * (std::fabs(v(0)) <= th);
							  } else
								  v(0) = 0. - emitter_length*(v(0) < 0);
						  }

						}

					  }
					}
				}
			break;
		  }

        default:
          {
            Assert(false, ExcNotImplemented());
          }
      }


	if (dim == 2)
	{
		GridGenerator::merge_triangulations( emitter_mesh, naca, triangulation);
		SetManifoldsAndBoundaries(triangulation);
	}
	else // caso 3D
	{
		Triangulation<2> aux;
		GridGenerator::merge_triangulations( emitter_mesh, naca, aux);
		const double width = 14.;
		unsigned int slices = 7;
		GridGenerator::extrude_triangulation(aux, slices, width, triangulation, true);
		SetManifoldsAndBoundaries(triangulation);
	}

// Switch per i vari tipi di emettitori
    switch (emitter_type)
      {
        case wire:
          {
        	const types::manifold_id wire_id = 1;

			if (dim==2)
			{
				const Point<2> wire_center(0.,0.);
				SphericalManifold<2> wire_manifold(wire_center);
				triangulation.set_manifold(wire_id, wire_manifold);
			}
			else
			{
	            Assert(false, ExcNotImplemented());
			}

			for (auto &face : triangulation.active_face_iterators()) {
			  if (face->at_boundary())
			  {
				  const Point<2> c = face->center();

				  if (c.square() < electrode_distance)
					  face->set_manifold_id(wire_id);
			  }
			}


			break;
          }

        case blade:
          {

  		    const double th = emitter_height*1e-6;
			const types::manifold_id	blade_id = 13;

			BladeGeometry<dim,2> blade;
			triangulation.set_manifold(blade_id, blade);

			for (auto &face : triangulation.active_face_iterators())
			{
			 if (face->at_boundary())
			  {
				  const Point<2> c = face->center();

				  if ( c.square() < electrode_distance ) {
					  Point<dim> &v = face->vertex(0);
					  if (std::fabs(c(0)-v(0)) > th && c(1) < - th)
						  face->set_manifold_id(blade_id);
				  }
			  }
			 }

            break;
          }

        case tip:
          {

  		    const double th = emitter_height*1e-6;
			const types::manifold_id	blade_up_id = 13;

			BladeGeometry<dim,1> blade_up;
			triangulation.set_manifold(blade_up_id, blade_up);

			const types::manifold_id	blade_down_id = 14;

			BladeGeometry<dim,-1> blade_down;
			triangulation.set_manifold(blade_down_id, blade_down);

			for (auto &face : triangulation.active_face_iterators())
			{
			 if (face->at_boundary())
			  {
				  const Point<2> c = face->center();

				  if ( c.square() < electrode_distance ) {
					  Point<dim> &v = face->vertex(0);
					  if (std::fabs(c(0)-v(0)) > th) {
						  if ( c(1) > th)
							  face->set_manifold_id(blade_up_id);
						  if (c(1) < -th)
							  face->set_manifold_id(blade_down_id);
					  }
				  }
			  }
			 }

            break;
          }

        case plate:
		  {
			break;
		  }

        default:
          {
            Assert(false, ExcNotImplemented());
          }
      }


}

// Definizione delle matrici e dei vettori per la soluzione
// Chiamata dopo ogni rifinimento della griglia per adeguare
// la dimensione
template <int dim>
void Problem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  constraints.clear(); // called multiple times
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close(); // rearranges constraints to optimize performance

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ true); // not overwriting matrix entries corresponding to constrained DoFs

  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
  laplace_matrix.reinit(sparsity_pattern);

  MatrixCreator::create_laplace_matrix(dof_handler,
                                       QGauss<dim>(fe.degree + 1),
                                       laplace_matrix);
}


// !!! Solver
template <int dim>
void Problem<dim>::solve()
{
 // Tolleranza e massimo numero di iterazioni
  unsigned int it_max = 1e+6;
  double tol = 1e-12;

  /*
   La tolleranza può anche essere definita in base alla norma l2
   del right hand side (che però nel nostro caso è assunto nullo):

   double tol = 1e-6* system_rhs.l2_norm();

   */

  // Tipo di solver: CG (Conjugate Gradient)
  SolverControl            solver_control(it_max, tol);
  SolverCG<Vector<double>> solver(solver_control);

  // Parametro di rilassamento e precondizionatore
  // Sono anche disponibili Jacobi e parametrizzazione LU (vedi sotto)
  double relaxation_parameter = 1.5;
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, relaxation_parameter);

  /* Alternative:
   PreconditionJacobi<SparseMatrix<double>> preconditioner;
   preconditioner.initialize(system_matrix);

   oppure:

   SparseILU<SparseMatrix<double>> preconditioner;
   preconditioner.initialize(system_matrix);
   */

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  // Stampa a schermo del numero di iterazioni richieste
  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;

  /* Il numero di iterazioni e il tempo CPU (sempre stampato a schermo)
   * sono utili per definire la bontà di un solutore, in particolare
   * studiando il loro andamento al variare degli N gradi di libertà per vedere
   * se seguono un andamento proporzionale a N o N^2, per esempio
   */

  constraints.distribute(solution); // computes the values of constrained nodes from the values of the unconstrained ones
}


// !!! Funzione per raffinamento griglia
template <int dim>
void Problem<dim>::refine_grid(const unsigned int min_grid_level,
        					   const unsigned int max_grid_level)
{
// Vettore errore (funziona solo con float, per quello che ho trovato)
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  // Stimatore dell'errore, studiato per l'equazione di Laplace
  KellyErrorEstimator<dim>::estimate(dof_handler,
                                     QGauss<dim - 1>(fe.degree + 1), // Formula di quadratura di Gauss di grado pari agli elementi finiti + 1
                                     {},
                                     solution,
                                     estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
												  std::min(0.1,0.01+0.025*cycle),  // % celle con errore maggiore, da rifinire
                                                  0.1); // % celle con errore minore, da riunire

  if (triangulation.n_levels() > max_grid_level)
    for (const auto &cell : triangulation.active_cell_iterators_on_level(max_grid_level))
      cell->clear_refine_flag();
  for (const auto &cell : triangulation.active_cell_iterators_on_level(min_grid_level))
    cell->clear_coarsen_flag();

  SolutionTransfer<dim> solution_trans(dof_handler);

   Vector<double> previous_solution;
   previous_solution = solution;
   triangulation.prepare_coarsening_and_refinement();
   solution_trans.prepare_for_coarsening_and_refinement(previous_solution);

   triangulation.execute_coarsening_and_refinement();
   setup_system();

   solution_trans.interpolate(previous_solution, solution);
   constraints.distribute(solution);
}



template <int dim>
void Problem<dim>::output_results()
{
	Point<2> sample;

	switch (emitter_type)
	{
		case wire:
		{
	        sample = {wire_radius,1e-3*emitter_height};
	        break;
		}

		case blade:
		{
	        sample = {1e-3*emitter_height,1e-3*emitter_height};
	        break;
		}

		case tip:
		{
	        sample = {0.,1e-3*emitter_height};
	        break;
		}

		case plate:
		{
	        sample = {1e-3*emitter_height,emitter_height/2};
	        break;
		}


		default:
		  {
			Assert(false, ExcNotImplemented());
		  }
	}

	/* Valore del potenziale in (1,0.2)
	 * (precedentemente usato come criterio di convergenza)
	 */
	float y = VectorTools::point_value(dof_handler, solution, {1.,0.2});
	std::cout << "   Potential at (1,0.2): " << y << std::endl;

	/* Se, per qualche motivo, potesse essere utile valutare la
	 convergenza sulla media del potenziale nel dominio,
	 questa è la funzione per approssimarla:

	 float x = VectorTools::compute_mean_value (dof_handler,QGauss<2>(fe.degree + 1),solution,0);

	 */

	// !!! Criterio di convergenza
	Tensor<1,dim>  E = VectorTools::point_gradient(dof_handler, solution, sample);
	float x = L2Norm(E);
	std::cout << "   Field at (" << sample[0] << "," << sample[1] << "): " << x << std::endl;

    if (cycle == 0) {
    	values.reinit(Nmax+2);
    	values(0) = 0;
    }

	values(cycle+1) = x;

	// Se la condizione è verificata, si considera si sia raggiunta convergenza
    if ( std::fabs(values(cycle) - x) <= conv_tol*std::fabs(x))
    	cycle = Nmax;




	// A convergenza raggiunta, stampa a schermo i risultati e viene scritta la soluzione su file
    if (cycle == Nmax) {

    	GradientPostprocessor<dim> gradient_postprocessor;
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution, "Potential");
        data_out.add_data_vector (solution, gradient_postprocessor);
        data_out.build_patches();

        std::ofstream output("solution-" + std::to_string(emitter_type) + ".vtu");
        data_out.write_vtu(output);

    	Tensor<1,dim>  E = VectorTools::point_gradient(dof_handler, solution, sample);
    	std::cout << "Electric field in ("<< sample[0] << "," << sample[1] << "): " << -E << ", magnitude: " << L2Norm(E) << std::endl;

    	E = VectorTools::point_gradient(dof_handler, solution, {sample[0],-sample[1]});
    	std::cout << "Electric field in ("<< sample[0] << "," << -sample[1] << "): " << -E << ", magnitude: " << L2Norm(E) << std::endl;

    	 E = VectorTools::point_gradient(dof_handler, solution, {electrode_distance,0.});
        std::cout << "Electric field before airfoil: " << -E << ", magnitude: " << L2Norm(E) << std::endl;

    	E = VectorTools::point_gradient(dof_handler, solution, {electrode_distance+airfoil_length,0.});
    	std::cout << "Electric field after airfoil: " << -E << ", magnitude: " << L2Norm(E) << std::endl;

        std::cout << std::endl;
    }
}



// Ciclo globale: definizione matrici e condizioni al contorno
template <int dim>
void Problem<dim>::run()
{
	const unsigned int max_refinement = 10;
	const unsigned int min_refinement = 0;

  while (cycle <= Nmax)
    {
      if (cycle == 0)
    	  create_mesh();
      else
        refine_grid(min_refinement, std::min(cycle + 1,max_refinement) );

  	  setup_system();

  	  // Forzante: imposta pari a 0
      VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), system_rhs);

      // Matrice: matrice di Laplace, moltiplicata per la permittività dell'aria
      system_matrix.copy_from(laplace_matrix);
      system_matrix *= 1.0006*eps0*1e-2;

      constraints.condense(system_matrix, system_rhs);

      // Condizioni al contorno
      {
        std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

        VectorTools::interpolate_boundary_values(dof_handler,
        											1, // Boundary corrispondente all'emettitore, definito sopra
													Functions::ConstantFunction<dim>(2.e+4), // Valore di potenziale all'emettitore (20 kV)
													emitter_boundary_values);

        VectorTools::interpolate_boundary_values(dof_handler,
        											2,  // Boundary corrispondente al collettore, definito sopra
													Functions::ZeroFunction<dim>(), // Valore di potenziale al collettore (0 V)
													collector_boundary_values);

        /* Le condizioni sopra sono condizioni di Dirichlet
        Su il restante bordo del dominio, dove non è imposta esplicitamente nessuna condizione,
        viene imposta automaticamente da deal.II una condizione di Neumann: gradiente nullo */

        MatrixTools::apply_boundary_values(emitter_boundary_values,
                                           system_matrix,
                                           solution,
  										   system_rhs);

        MatrixTools::apply_boundary_values(collector_boundary_values,
                                           system_matrix,
                                           solution,
  										   system_rhs);
      }


      std::cout << std::endl
      	  	  	<< "Cycle " << cycle << ':' << std::endl
				<< "   Number of active cells:       " << triangulation.n_active_cells() << std::endl
      	  	    << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

      solve();
      output_results();

	  std::cout << "   Elapsed CPU time: " << timer.cpu_time() << " seconds.\n";

      cycle++;
    }
}


/* Ciclo main per i vari emettitori
 * (generalmente commento o "s-commento"
 * a seconda di quale emettitore voglio usare:
 * vengono simulati tutti in sequenza)
 */
int main()
{
  try
    {
	  /* Il numero in <> indica la dimensione del problema
	   la parola dopo i due punti il tipo di emettitore  */

      Problem<2> wire_poisson_2d(Problem<2>::wire); // Definisco il tipo di problema...
      wire_poisson_2d.run(); // ... e lo simulo

	  Problem<2> tip_poisson_2d(Problem<2>::tip);
	  tip_poisson_2d.run();

	  Problem<2> blade_poisson_2d(Problem<2>::blade);
      blade_poisson_2d.run();

      Problem<2> plate_poisson_2d(Problem<2>::plate);
      plate_poisson_2d.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1; // Report an error
    }
  catch (...)
    {
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
