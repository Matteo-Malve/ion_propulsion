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

  

  // ------------------------------------------------------------      
  // TABLE
  // ------------------------------------------------------------      
  
  convergence_table.add_value("cycle", cycle);
  convergence_table.add_value("cells", triangulation.n_active_cells());
  Evaluation::PointValueEvaluation<dim> postprocessor(EVALUATION_POINT);
  double computed_value = postprocessor(primal_dof_handler,uh);
  cout<<"      Point Value:   "<< computed_value << endl;
  convergence_table.add_value("Point value", computed_value);
 

}


// #######################################
// Template initialization
// #######################################
template class Problem<2>;