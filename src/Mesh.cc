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

  Problem();

  void run();

private:
  void create_mesh(const double mesh_height, const double electrode_distance,
		           const double wire_radius, const double collector_height);
  void setup_system();
  void solve();
  void refine_grid(const unsigned int min_grid_level, // minimo livello di raffinamento griglia
          	  	  const unsigned int max_grid_level); // massimo livello di raffinemento griglia
  void output_results(const double wire_radius);

  Triangulation<dim> triangulation;

  FE_Q<dim>       fe;
  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints;

  SparseMatrix<double> system_matrix;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> laplace_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;

  Timer timer;

  unsigned int cycle = 0;
  const unsigned int Nmax = 10; // massimo numero di cicli

  Vector<float> values;
  const float conv_tol = 1e-4; // tolleranza globale
};


// Funzione per stampare a schermo alcune caratteristiche della griglia
template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename)
{
  std::cout << "Mesh info:" << std::endl // Comando per stampare a schermo, utilissimo per il debugging
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;

    std::cout << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        std::cout << pair.first << " (" << pair.second << " times) ";
      }
    std::cout << std::endl;

    std::map<types::manifold_id, unsigned int> manifold_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        manifold_count[face->manifold_id()]++;

    std::cout << " manifold indicators: ";
    for (const std::pair<const types::manifold_id, unsigned int> &pair :
         manifold_count)
      {
        std::cout << pair.first << " (" << pair.second << " times) ";
      }
    std::cout << std::endl;
  }
  std::ofstream out(filename);
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << " written to " << filename << std::endl << std::endl;
}


// Funzione per calcolo norma di tensore
template <int dim>
double L2Norm(const Tensor<1,dim> &input)
{
	double magnitude = 0.;
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


//  Funzione per creare la griglia
template <int dim>
void CreateGrid( Triangulation<dim> &mesh, const double mesh_height,
				 const double electrode_distance, const double wire_radius)
{
	Assert(dim == 2, ExcNotImplemented()); // Serve solo per dare errore se si prova con dim = 3

	const Point<2> top_left(-electrode_distance-wire_radius, mesh_height); // Coordinate angolo in alto a sinistra
	const Point<2> bottom_right(electrode_distance+wire_radius, 0.); // Coordinate angolo in basso a destra

	unsigned int a = (int) std::ceil( (bottom_right(0) - top_left(0))/wire_radius ); // Numero di suddivisioni lungo x
	unsigned int b = (int) std::ceil( (top_left(1) - bottom_right(1))/wire_radius/2 ); // Numero di suddivisioni lungo y

	// Griglia rettangolare:
	GridGenerator::subdivided_hyper_rectangle( mesh, {a,b}, top_left, bottom_right);


	// Ciclo su tutte le facce per spostare quelle dove si trova l'elettrodo
	  for (auto &face : mesh.active_face_iterators())
	  {
		  if (face->at_boundary()) // Se la faccia è sul bordo ...
		  {
			  const Point<2> c = face->center();
			  if ( std::fabs(c[0]) <= wire_radius && c[1] < wire_radius*1e-3 ) // ... e il centro dista da 0,0 meno del raggio dell'emettitore
			  {
				for (const auto i : face->vertex_indices()) // Ciclo sui vertici della cella
				{
				  Point<2> &v = face->vertex(i);
				  v(1) = std::max(0. , std::sqrt(pow(wire_radius,2) - pow(v(0),2) ) ); // impongo y = sqrt(r^2 - x^2)
				}
		  }
		  }
	 }


}


// Funzione per assegnare i vari manifold e il numero identificativo dei boundaries
template <int dim>
void SetManifoldsAndBoundaries(Triangulation<dim> &mesh, const double collector_height,
								const double electrode_distance, const double wire_radius)
{
	Assert(dim == 2, ExcNotImplemented());

	// Boundaries
	const types::boundary_id collector_id = 2;
	const types::boundary_id emitter_id = 1;

	// Manifolds:
	 const types::manifold_id     wire_id = 1;

	 const Point<2>       center(0., 0.);
	 SphericalManifold<dim> circle(center);
	 mesh.set_manifold(wire_id, circle);

	  for (auto &face : mesh.active_face_iterators())
	  {
		  if (face->at_boundary())
		  {
			  const Point<dim> c = face->center();

			  if ( (c[1] > 0) && (c[1] <= collector_height) && (c[0] > electrode_distance) )
				  face->set_boundary_id(collector_id);
			  else if ( (c[1] > 0) && (std::fabs(c[0]) < 2*wire_radius) && (c[1] <= wire_radius*2)) {
				  face->set_boundary_id(emitter_id);
				  face->set_manifold_id(wire_id);
				  //std::cout << " Set wire manifold in: " << c << std::endl;
			  }
		  }
	  }
}

// !!! Scelta degli elementi finiti: lineari o quadratici
template <int dim>
Problem<dim>::Problem()
  : fe(1) // elemeenti lineari (1) o quadratici (2)
  , dof_handler(triangulation)
{}


/* Creazione della griglia computazionale
Tutorial griglie: step-49 su deal.ii */
template <int dim>
void Problem<dim>::create_mesh(const double mesh_height, const double electrode_distance,
								const double wire_radius, const double collector_height)
{
	CreateGrid(triangulation, mesh_height, electrode_distance, wire_radius);

	SetManifoldsAndBoundaries(triangulation, std::min(mesh_height,collector_height), electrode_distance, wire_radius);

	print_mesh_info(triangulation, "grid.svg");

	triangulation.refine_global(1);

	print_mesh_info(triangulation, "refined_grid.svg");

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

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ true);

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
  unsigned int it_max = 1000;
  double tol = 1e-14;

  /*
   La tolleranza può anche essere definita in base alla norma l2
   del right hand side (che però nel nostro caso è assunto nullo):

   double tol = 1e-6 * system_rhs.l2_norm();

   */

  // Tipo di solver: CG (Conjugate Gradient)
  SolverControl            solver_control(it_max, tol);
  SolverCG<Vector<double>> solver(solver_control);

  // Parametro di rilassamento e precondizionatore
  // Sono anche disponibili Jacobi e parametrizzazione LU (vedi sotto)
  double relaxation_parameter = 1.2; // Questo parametro può variare nell'intervallo [1,2)
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
   sono utili per definire la bontà di un solutore, in particolare
   studiando il loro andamento al variare degli N gradi di libertà
   della griglia per vedere se seguono un andamento proporzionale
   a N o N^2, per esempio
   */

  constraints.distribute(solution);
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
												  0.05 + 0.02*cycle,  // % celle con errore maggiore (da rifinire)
                                                  0.1); // % celle con errore minore (da "ri-unire")

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
void Problem<dim>::output_results(const double wire_radius)
{
	// Inizializzazione:
    if (cycle == 0) {
    	values.reinit(Nmax+2);
    	values(0) = 0;
    }

	/* Valore del potenziale in (1,0.2)
	 (precedentemente usato come criterio di convergenza)

	double y = VectorTools::point_value(dof_handler, solution, {1.,0.2});
	std::cout << "   Potential at (1,0.2): " << y << std::endl;

    */

	/* Se, per qualche motivo, potesse essere utile valutare la
	 convergenza sulla media del potenziale nel dominio,
	 questa è la funzione per approssimarla:

	 double x = VectorTools::compute_mean_value (dof_handler,QGauss<2>(fe.degree + 1),solution,0);

	 */

	// !!! Criterio di convergenza
	Point<dim> sample(1.,0.2);

	Tensor<1,dim>  E = VectorTools::point_gradient(dof_handler, solution, sample);
	double x = L2Norm(E);

	std::cout << "   Field magnitude at (" << sample[0] << "," << sample[1] << "): " << x << std::endl;

	values(cycle+1) = x;

	// Se la condizione è verificata, si considera si sia raggiunta convergenza
    if ( std::fabs(values(cycle) - x) <= conv_tol*std::fabs(x))
    	cycle = Nmax;

    // Scrittura su file della soluzione ad ogni ciclo:
	GradientPostprocessor<dim> gradient_postprocessor;
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "Potential");
    data_out.add_data_vector (solution, gradient_postprocessor);
    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
    data_out.write_vtu(output);


	// A convergenza raggiunta, stampa a schermo i risultati
    if (cycle == Nmax) {

    	Point<dim> sample(wire_radius, 0.);
    	Tensor<1,dim>  E = VectorTools::point_gradient(dof_handler, solution, sample);
    	std::cout << "Electric field in ("<< sample[0] << "," << sample[1] << "): " << -E << ", magnitude: " << L2Norm(E) << std::endl;

        std::cout << std::endl;
    }
}



// Ciclo globale: definizione matrici e condizioni al contorno
template <int dim>
void Problem<dim>::run()
{
	// Variabili mesh
	const double mesh_height = 4.; // mm
	const double electrode_distance = 2.; // mm
	const double wire_radius = 0.025; // mm
	const double collector_height = 1.2; // mm

	// Raffinamento griglia massimo e minimo
	const unsigned int max_refinement = 20;
	const unsigned int min_refinement = 0;

	std::cout << "Simulating for WIRE of radius: " << wire_radius
			  << " mm in " << dim << "D" << std::endl
			  << "Mesh: " << 2*electrode_distance << " mm x " << mesh_height << " mm"
			  << std::endl;

  while (cycle <= Nmax)
    {

      if (cycle == 0)
    	  create_mesh(mesh_height, electrode_distance, wire_radius, collector_height);
      else
        refine_grid(min_refinement, std::min(cycle + 1,max_refinement) );

  	  setup_system();

  	  // Forzante: imposta pari a 0
      VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), system_rhs);

      // Permittività elettrica del vuoto:
      const double eps0 = 8.854*1e-12; // [F/m]

      // Matrice: matrice di Laplace, moltiplicata per la permittività dell'aria
      system_matrix.copy_from(laplace_matrix);
      system_matrix *= 1.0006*eps0*1e-3;

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
      output_results(wire_radius);

	  std::cout << "   Elapsed CPU time: " << timer.cpu_time() << " seconds.\n";

      cycle++;
    }
}


// Ciclo main
int main()
{
  try
    {
	  // Il numero in <> indica la dimensione del problema

      Problem<2> wire_poisson_2d; // Definisco il tipo di problema...
      wire_poisson_2d.run(); // ... e lo simulo

      // Non preoccupatevi del caso 3D: non funziona...
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
