#include "problem.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

using namespace dealii;
using std::cout;
using std::endl;

template <int dim>
double Problem<dim>::compute_averaged_error() const {
    // Initialize error accumulator and point counter
    double total_error = 0.0;
    unsigned int total_points = 0;

    // Select a quadrature rule (e.g., QGauss) for integration
    const unsigned int quadrature_order = 2;  // Adjust this based on your needs
    dealii::QGauss<dim> quadrature_formula(quadrature_order);

    // Initialize FEValues to evaluate the solution at quadrature points
    dealii::FEValues<dim> fe_values(primal_fe, quadrature_formula,
                                    dealii::update_values | dealii::update_quadrature_points);

    // Iterate over each cell
    for (const auto &cell : primal_dof_handler.active_cell_iterators()) {
        // Reinitialize FEValues for the current cell
        fe_values.reinit(cell);

        // Extract solution values at quadrature points on the current cell
        std::vector<double> uh_values(fe_values.n_quadrature_points);
        fe_values.get_function_values(uh, uh_values);

        // Collect quadrature points and evaluate exact solution at all points in one call
        std::vector<dealii::Point<dim>> quadrature_points = fe_values.get_quadrature_points();
        std::vector<double> uex_values(fe_values.n_quadrature_points);
        exact_solution_function.value_list(quadrature_points, uex_values);

        // Loop over quadrature points on this cell to calculate errors
        for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
            // Compute error at this quadrature point
            double error = std::fabs(uh_values[q] - uex_values[q]);

            // Accumulate the error and increase the point counter
            total_error += error;
            ++total_points;
        }
    }

    // Compute the average error by dividing the accumulated error by the number of points
    double averaged_error = total_error / total_points;
    return averaged_error;
}


template <int dim>
double Problem<dim>::localized_average_error(Point<dim> center_point, double radius) const {
    // Initialize error accumulator and point counter
    double total_error = 5.0;
    unsigned int total_points = 0;
    //cout<<"      radius = "<<radius<<endl;
    //cout<<"      center point = "<<center_point(0)<<" "<<center_point(1)<<endl;

    // Select a quadrature rule (e.g., QGauss) for integration
    const unsigned int quadrature_order = 2;  // Adjust this based on your needs
    QGauss<dim> quadrature_formula(quadrature_order);

    // Initialize FEValues to evaluate the solution at quadrature points
    FEValues<dim> fe_values(primal_fe, quadrature_formula,
                                    update_values | update_quadrature_points);
    
    // Iterate over each cell
    for (const auto &cell : primal_dof_handler.active_cell_iterators()) {
        // Check if the cell center is within the radius of the specified center point
        Point<dim> cell_center = cell->center();
        double distance = cell_center.distance(center_point);
        if(distance<5.*radius)
          //cout<<"      distance: "<<distance<<endl;
        if (distance <= radius) {
            // Reinitialize FEValues for the current cell
            fe_values.reinit(cell);

            // Extract solution values at quadrature points on the current cell
            std::vector<double> uh_values(fe_values.n_quadrature_points);
            fe_values.get_function_values(uh, uh_values);

            // Loop over quadrature points on this cell
            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
                // Exact solution at this quadrature point
                double uex_value = exact_solution_function.value(fe_values.quadrature_point(q));

                // Compute error at this quadrature point
                double error = std::fabs(uh_values[q] - uex_value);

                // Accumulate the error and increase the point counter
                total_error += error;
                ++total_points;
            }
        }
    }
    // Compute the average error by dividing the accumulated error by the number of points
    double averaged_error = (total_points > 0) ? (total_error / total_points) : 0.0;
    return averaged_error;
}

template <int dim>
double Problem<dim>::compute_localized_L2_error(const Point<dim> &center_point, const double radius) {
    // Set up an error vector
    Vector<double> difference_per_cell(triangulation.n_active_cells());

    // Compute the difference between the finite element solution and the exact solution
    const QGauss<dim> quadrature_formula(2 * primal_dof_handler.get_fe().degree + 1);
    VectorTools::integrate_difference(
        primal_dof_handler,
        uh,
        exact_solution_function,
        difference_per_cell,
        quadrature_formula,
        VectorTools::L2_norm);

    // Loop through cells and filter the ones within the radius
    double localized_L2_error = 0.0;
    for (const auto &cell : triangulation.active_cell_iterators()) {
        const Point<dim> cell_center = cell->center();
        if (cell_center.distance(center_point) < radius) {
            localized_L2_error += difference_per_cell[cell->active_cell_index()];
        }
    }

    return std::sqrt(localized_L2_error);
}

template <int dim>
double Problem<dim>::compute_localized_H1_error(const Point<dim> &center_point, const double radius) {
    // Set up an error vector
    Vector<double> difference_per_cell(triangulation.n_active_cells());

    // Compute the difference between the finite element solution and the exact solution
    const QGauss<dim> quadrature_formula(7);
    VectorTools::integrate_difference(primal_dof_handler,
                                      uh,
                                      exact_solution_function,
                                      difference_per_cell,
                                      quadrature_formula,
                                      VectorTools::H1_norm);

    // Loop through cells and filter the ones within the radius
    double localized_H1_error_squared = 0.0;
    for (const auto &cell : triangulation.active_cell_iterators()) {
        const Point<dim> cell_center = cell->center();
        if (cell_center.distance(center_point) < radius) {
            localized_H1_error_squared += difference_per_cell[cell->active_cell_index()];
        }
    }

    // Return the square root of the accumulated error
    return std::sqrt(localized_H1_error_squared);
}

template <int dim>
double Problem<dim>::compute_L2_error() {

    // Set up an error vector
    Vector<double> difference_per_cell(triangulation.n_active_cells());

    // Compute the difference between the finite element solution and the exact solution
    const QGauss<dim> quadrature_formula(2*primal_dof_handler.get_fe().degree + 1);
    VectorTools::integrate_difference(
        primal_dof_handler,
        uh,
        exact_solution_function,
        difference_per_cell,
        quadrature_formula,
        VectorTools::L2_norm);

    const double L2_error = VectorTools::compute_global_error(triangulation,
                                                              difference_per_cell,
                                                              VectorTools::L2_norm);
    return L2_error;
}

template <int dim>
double Problem<dim>::compute_H1_error() {

    // Set up an error vector
    Vector<double> difference_per_cell(triangulation.n_active_cells());

    // Compute the difference between the finite element solution and the exact solution
    const QGauss<dim> quadrature_formula(7);
    VectorTools::integrate_difference(primal_dof_handler,                       // Different function for other Mappings
                                      uh,
                                      exact_solution_function,
                                      difference_per_cell,
                                      quadrature_formula,
                                      VectorTools::H1_norm);

    const double H1_error = VectorTools::compute_global_error(triangulation,
                                                              difference_per_cell,
                                                              VectorTools::H1_norm);  // Different function for other Mappings
    return H1_error;
}

template <int dim>
void Problem<dim>::test_convergence(){

  cout<<"   Convergence test:"<<endl;

  
  // ------------------------------------------------------------      
  // Average on all mesh points
  // ------------------------------------------------------------      

  /*double average_error = compute_averaged_error();
  cout<<"      Average error:           "<< average_error << endl;
  average_errors.push_back(average_error);*/

  // ------------------------------------------------------------      
  // L2 error
  // ------------------------------------------------------------      
  double L2_error = 0.0;
  //L2_error = compute_L2_error();
  cout<<"      L2 error:                "<< L2_error << endl;

  // ------------------------------------------------------------      
  // H1 error
  // ------------------------------------------------------------      
  double H1_error = 0.0;
  //H1_error = compute_H1_error();
  cout<<"      H1 error:                "<< H1_error << endl;

	// ------------------------------------------------------------      
  // localized L2 error
  // ------------------------------------------------------------      

  //double loc_L2_error = compute_localized_L2_error(EVALUATION_POINT, EVALUATION_RADIUS);
  //cout<<"      Localized L2 error:      "<< loc_L2_error << endl;

  // ------------------------------------------------------------      
  // localized H1 error
  // ------------------------------------------------------------      

  //double loc_H1_error = compute_localized_H1_error(EVALUATION_POINT, EVALUATION_RADIUS);
	//cout<<"      Localized H1 error:      "<< loc_H1_error << endl;

  // ------------------------------------------------------------      
  // Localized average in a ball
  // ------------------------------------------------------------      
  double loc_average_error = 0.0;
  if(GOAL_FUNCTIONAL == "AreaEvaluation")
  {
    loc_average_error = localized_average_error(EVALUATION_POINT, EVALUATION_RADIUS);
    cout<<"      Localized average error: "<< loc_average_error << endl;
  }
  

  // ------------------------------------------------------------      
  // Evaluation in target point
  // ------------------------------------------------------------     
  double error_target_point = 0.;
  if(GOAL_FUNCTIONAL == "PointValue")
  {
    //double exact_value = exact_solution_function.value(EVALUATION_POINT);
    double exact_value = 1.58111e-07;
    Evaluation::PointValueEvaluation<dim> postprocessor(EVALUATION_POINT);
    double computed_value = postprocessor(primal_dof_handler,uh);
    error_target_point = std::fabs(exact_value-computed_value);
    cout<<"      Exact value:             "<< exact_value << endl;
    cout<<"      Computed value:          "<< computed_value << endl;
    cout<<"      Error at target point:   "<< error_target_point << endl;
  }
	// ------------------------------------------------------------      
  // Evaluation of d(phi)/d(y) at target point
  // ------------------------------------------------------------      
  double error_dy_target_point = 0.;
  if(GOAL_FUNCTIONAL == "PointYDerivative")
  {
    double exact_value = exact_solution_function.gradient(EVALUATION_POINT)[1];
    Evaluation::PointYDerivativeEvaluation<dim> postprocessor(EVALUATION_POINT);
    double computed_value = postprocessor(primal_dof_handler,uh);
    error_dy_target_point = std::fabs(exact_value-computed_value);
    cout<<"      d(phi)/d(y) exact:      "<< exact_value << endl;
    cout<<"      d(phi)/d(y) computed:   "<< computed_value << endl;
    cout<<"      Error of d(phi)/d(y) at target point:   "<< error_dy_target_point << endl;
  }

	// ------------------------------------------------------------      
  // Emitter Flux evaluation
  // ------------------------------------------------------------ 
	double flux_error = 0;
	/*if(GOAL_FUNCTIONAL == "BoundaryFluxEvaluation"){
		double exact_flux = exact_solution_function.emitter_flux();
		cout<<"      Flux exact:             "<< exact_flux <<endl;
		Evaluation::FluxEvaluation<dim> postprocessor;
		double computed_flux = postprocessor(primal_dof_handler, uh);
		cout<<"      Flux computed:          "<< computed_flux <<endl;
		flux_error = std::fabs(exact_flux-computed_flux);
		cout<<"      Flux error:             "<< flux_error <<endl;
	}*/
  // ------------------------------------------------------------      
  // TEXT OUTPUT
  // ------------------------------------------------------------      

  // Prepare CSV file for writing
  std::ofstream csv_file;
  std::string file_name = TEST_NAME + "-convergence_data.csv";

  if(cycle==0){
    csv_file.open(file_name);

    // Write the header
    csv_file << "cycle,num_cells,DoFs,estimated_error,exact_error\n";

    // Close the file
    csv_file.close();

  } else{

    // Open the file in append mode to add new lines
    csv_file.open(file_name, std::ios::app);

    // Write the current cycle's data
    csv_file << cycle << ","    
             << triangulation.n_active_cells() << ","
             << primal_dof_handler.n_dofs() << ","
             << estimated_error << ",";
    
    if(GOAL_FUNCTIONAL == "PointValue")
      csv_file << error_target_point << "\n";
    else if(GOAL_FUNCTIONAL == "PointYDerivative")
      csv_file << error_dy_target_point << "\n";
    else if(GOAL_FUNCTIONAL == "AreaEvaluation")
      csv_file << loc_average_error << "\n";
    else if(GOAL_FUNCTIONAL == "BoundaryFluxEvaluation")
      csv_file << flux_error << "\n";
    else{
      cout<< "ERROR: Goal Functional strategy required is unkown."<<endl;
      abort();
    }

    // Close the file
    csv_file.close();
  }     
  

  // ------------------------------------------------------------      
  // TABLE
  // ------------------------------------------------------------      
  
  convergence_table.add_value("cycle", cycle);
  convergence_table.add_value("cells", triangulation.n_active_cells());
  convergence_table.add_value("DoFs", primal_dof_handler.n_dofs());
  convergence_table.add_value("estimated_error", estimated_error);

  if(GOAL_FUNCTIONAL == "PointValue")
    convergence_table.add_value("ex. point err", error_target_point);
  else if(GOAL_FUNCTIONAL == "PointYDerivative")
    convergence_table.add_value("ex. pointDy err", error_dy_target_point);
  else if(GOAL_FUNCTIONAL == "AreaEvaluation")
    convergence_table.add_value("ex. point err", loc_average_error);
  else if(GOAL_FUNCTIONAL == "BoundaryFluxEvaluation")
    convergence_table.add_value("ex. flux err", flux_error);
  else{
    cout<< "ERROR: Goal Functional strategy required is unkown."<<endl;
    abort();
  }

  convergence_table.add_value("L2", L2_error);
  convergence_table.add_value("H1", H1_error);
  
}


// #######################################
// Template initialization
// #######################################
template class Problem<2>;