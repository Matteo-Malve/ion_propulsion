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


/* Helper function, gives errors on less recent versions of deal.ii
 * It outputs a .vtu with a field dedicated to the boundary ids.
 * It is helpful to visualize if you've correctly set the desired BCs.

template <int dim>
void check_boundary_ids(const Triangulation<dim> &triangulation) {
    cout<<"Starting check on Boundary ids..."<<endl;
    DataPostprocessors::BoundaryIds <dim> boundary_ids;
    DataOutFaces<dim> data_out_faces; // requires: #include <deal.II/numerics/data_out_faces.h>
    FE_Q <dim> dummy_fe(1);

    DoFHandler <dim> dummy_dof_handler(triangulation);
    dummy_dof_handler.distribute_dofs(dummy_fe);

    Vector<double> dummy_solution(dummy_dof_handler.n_dofs());

    data_out_faces.attach_dof_handler(dummy_dof_handler);
    data_out_faces.add_data_vector(dummy_solution, boundary_ids);
    data_out_faces.build_patches();

    std::ofstream out("boundary_ids.vtu");
    data_out_faces.write_vtu(out);
    cout<<"   Boundary ids written to boundary_ids.vtu"<<endl
        <<"   Check the file\n\n";     
}  */




// #######################################
// Template initialization
// #######################################

template double L2Norm(const Tensor<1,2> &input);
template void ionization_area(const Triangulation<2> &triangulation, const DoFHandler<2> &dof_handler, const Vector<double> &solution);
//template void check_boundary_ids(const Triangulation<2> &triangulation);