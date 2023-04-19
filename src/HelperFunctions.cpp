#include "../include/HelperFunctions.h"
//#include "HelperFunctions.h"

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

template <int dim>
double L2Norm(const Tensor<1,dim> &input)
{
    double magnitude = 0.;
    for (unsigned int i=0; i<dim; ++i)
        magnitude += input[i]*input[i];

    return std::sqrt(magnitude);
}


// #######################################
// Template initialization
// #######################################
template void print_mesh_info(const Triangulation<2> &triangulation,
                              const std::string &       filename);
template double L2Norm(const Tensor<1,2> &input);
