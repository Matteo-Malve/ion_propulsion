#include "HelperFunctions.h"
#include "../Evaluation.h"

template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename)
{
    std::cout << "\nMesh info:" << std::endl
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
    //const double E_soglia = datafile("Ionization/E_soglia",4e+15);
    const double E_soglia = 2e+3;

    for (const auto &cell : triangulation.active_cell_iterators()){
        Point<dim> c = cell->center();
        //Evaluation::PointValueEvaluation<dim> postprocessor(c);
        //double x_ = postprocessor(dof_handler, solution);
        Tensor<1,dim>  E_ = VectorTools::point_gradient(dof_handler, solution, c);
        double x = L2Norm(E_);
        if (x > E_soglia)
            cell->set_material_id(2);
    }    

    std::ofstream out("ionization_area.vtu");
    GridOut       grid_out;
    grid_out.write_vtu(triangulation, out);
}

/* template <int dim>
void ckeck_boundary_ids(const Triangulation<dim> &triangulation) {
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
template void print_mesh_info(const Triangulation<2> &triangulation,
                              const std::string &       filename);
template double L2Norm(const Tensor<1,2> &input);
template void ionization_area(const Triangulation<2> &triangulation, const DoFHandler<2> &dof_handler, const Vector<double> &solution);
//template void ckeck_boundary_ids(const Triangulation<2> &triangulation);