#include "problem.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

using namespace dealii;
using std::cout;
using std::endl;
namespace plt = matplotlibcpp;

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
    cout<<"      radius = "<<radius<<endl;
    cout<<"      center point = "<<center_point(0)<<" "<<center_point(1)<<endl;

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
void Problem<dim>::test_convergence(){
  Point<dim> sensor_1(0.0002625, 0.0005);
  Point<dim> sensor_2(-0.00025, 0.0005);
  Point<dim> sensor_3(0.0031, 0.0035);
  Point<dim> sensor_4(-0.0004, -0.0013);

  cout<<"   Convergence test:"<<endl;
  double uh_at_sensor_1 = VectorTools::point_value(primal_dof_handler, uh, sensor_1);
  double uex_at_sensor_1 = exact_solution_function.value(sensor_1);
  double abs_err_1 = std::fabs(uh_at_sensor_1-uex_at_sensor_1);
  /*std::cout << "      Sensor 1:" << endl
            << "         uh      =  " << uh_at_sensor_1 << endl
            << "         u_ex    =  " << uex_at_sensor_1 << endl
            << "         abs_err =  " << abs_err_1 <<endl;*/
  
  double uh_at_sensor_2 = VectorTools::point_value(primal_dof_handler, uh, sensor_2);
  double uex_at_sensor_2 = exact_solution_function.value(sensor_2);
  double abs_err_2 = std::fabs(uh_at_sensor_2-uex_at_sensor_2);
  /*std::cout << "      Sensor 2:" << endl
            << "         uh      =  " << uh_at_sensor_2 << endl
            << "         u_ex    =  " << uex_at_sensor_2 << endl
            << "         abs_err =  " << abs_err_2 <<endl;*/
  
  double uh_at_sensor_3 = VectorTools::point_value(primal_dof_handler, uh, sensor_3);
  double uex_at_sensor_3 = exact_solution_function.value(sensor_3);
  double abs_err_3 = std::fabs(uh_at_sensor_3-uex_at_sensor_3);
  /*std::cout << "      Sensor 3:" << endl
            << "         uh      =  " << uh_at_sensor_3 << endl
            << "         u_ex    =  " << uex_at_sensor_3 << endl
            << "         abs_err =  " << abs_err_3 <<endl;*/
  
  double uh_at_sensor_4 = VectorTools::point_value(primal_dof_handler, uh, sensor_4);
  double uex_at_sensor_4 = exact_solution_function.value(sensor_4);
  double abs_err_4 = std::fabs(uh_at_sensor_4-uex_at_sensor_4);
  /*std::cout << "      Sensor 4:" << endl
            << "         uh      =  " << uh_at_sensor_4 << endl
            << "         u_ex    =  " << uex_at_sensor_4 << endl
            << "         abs_err =  " << abs_err_4 <<endl;*/
  
  errors_sensor_1.push_back(abs_err_1);
  errors_sensor_2.push_back(abs_err_2);
  errors_sensor_3.push_back(abs_err_3);
  errors_sensor_4.push_back(abs_err_4);
  double average_error = compute_averaged_error();
  cout<<"      Average error: "<< average_error << endl;
  average_errors.push_back(average_error);
  evaluation_point = Point<dim>(0.00025, 0.0005);  
  double loc_average_error = localized_average_error(evaluation_point, R/3.);
  cout<<"      Localized average error: "<< loc_average_error << endl;
  localized_average_errors.push_back(loc_average_error);

  // Prepare CSV file for writing
  std::ofstream csv_file;
  std::string file_name = TEST_NAME + "-convergence_data.csv";
  csv_file.open(file_name);
  
  // Write the header
  csv_file << "num_cells,errors_sensor_1,errors_sensor_2,errors_sensor_3,errors_sensor_4,average_errors,localized_average_errors\n";
  
  // Write each cycle's data
  for (size_t i = 0; i < num_cells.size(); ++i) {
      csv_file << num_cells[i] << ","
               << errors_sensor_1[i] << ","
               << errors_sensor_2[i] << ","
               << errors_sensor_3[i] << ","
               << errors_sensor_4[i] << ","
               << average_errors[i] << ","
               << localized_average_errors[i] << "\n";
  }
  
  // Close the CSV file
  csv_file.close();

  // Plot data
  plt::figure_size(800, 600);
  plt::clf();
  
  /*plt::named_plot("Sensor 1", cycles, errors_sensor_1, "r-o");
  plt::named_plot("Sensor 2", cycles, errors_sensor_2, "g-o");
  plt::named_plot("Sensor 3", cycles, errors_sensor_3, "b-o");
  plt::named_plot("Sensor 4", cycles, errors_sensor_4, "k-o");
  plt::named_plot("average", cycles, average_errors, "y-o");
  plt::named_plot("localized average", cycles, localized_average_errors, "y:o");*/

  plt::named_loglog("Sensor 1", num_cells, errors_sensor_1, "r-o");
  plt::named_loglog("Sensor 2", num_cells, errors_sensor_2, "g-o");
  plt::named_loglog("Sensor 3", num_cells, errors_sensor_3, "b-o");
  plt::named_loglog("Sensor 4", num_cells, errors_sensor_4, "k-o");
  plt::named_loglog("average", num_cells, average_errors, "y-o");
  plt::named_loglog("localized average", num_cells, localized_average_errors, "y:o");


  //plt::xlabel("Cycle");
  plt::xlabel("Number of Cells");
  plt::ylabel("Error");
  plt::title("Convergence error");
  plt::legend();
  plt::grid(true);

  // Save the plot
  plt::save(TEST_NAME+"-convergence_test.png");
}


// #######################################
// Template initialization
// #######################################
template class Problem<2>;