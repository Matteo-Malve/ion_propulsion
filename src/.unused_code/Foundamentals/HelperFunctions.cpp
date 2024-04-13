#include "HelperFunctions.h"
#include "../Evaluation.h"

template <int dim>
double L2Norm(const Tensor<1,dim> &input)
{
    double magnitude = 0.;
    for (unsigned int i=0; i<dim; ++i)
        magnitude += input[i]*input[i];
    return std::sqrt(magnitude);
}

template <int dim>
void ionization_area(const Triangulation<dim> &triangulation, const DoFHandler<dim> &dof_handler, const Vector<double> &solution) {
    const double E_threshold = redefined_4_datafile("electrical_field_intensity_threshold",2e+5);
    for (const auto &cell : triangulation.active_cell_iterators()){
        Point<dim> c = cell->center();
        Tensor<1,dim>  E_ = VectorTools::point_gradient(dof_handler, solution, c);
        double x = L2Norm(E_);
        if (x > E_threshold)
            cell->set_material_id(2);
    }    
    // At each refinement loop it will overwrite itself. Only last one remains.
    std::ofstream out("ionization_area.vtu");
    GridOut       grid_out;
    grid_out.write_vtu(triangulation, out);
}

// #######################################
// Template initialization
// #######################################

template double L2Norm(const Tensor<1,2> &input);
template void ionization_area(const Triangulation<2> &triangulation, const DoFHandler<2> &dof_handler, const Vector<double> &solution);
