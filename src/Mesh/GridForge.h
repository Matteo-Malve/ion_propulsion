#ifndef ION_PROPULSION_GRIDFORGE_H
#define ION_PROPULSION_GRIDFORGE_H

#include "../includes&parameters_setup.h"



template <int dim>
void CreateGrid(Triangulation<dim> &mesh);

template <int dim>
void LoadSecondGrid(Triangulation<dim> &mesh);

template <int dim>
void CreateInitialGrid( Triangulation<dim> &mesh, const double mesh_height,
                 const double electrode_distance, const double wire_radius);

template<int dim>
void FirstRefineGrid(Triangulation<dim> &);

template <int dim>
void CurveCells( Triangulation<dim> &mesh, const double mesh_height,
                        const double electrode_distance, const double wire_radius);

/* ----- Useless now that we fix boundary ids in gmsh
template <int dim>
void SetManifoldsAndBoundaries(Triangulation<dim> &mesh, const double collector_height,
                               const double electrode_distance, const double wire_radius);
*/
/*
template<int dim>
void CreateSecondGrid(Triangulation<dim> &);

template<>
void CreateSecondGrid(Triangulation<2> & triangulation){
        const Point<2> center(1, 0);
        const double inner_radius = 0.5, outer_radius = 1.0;
        GridGenerator::hyper_shell(
                triangulation, center, inner_radius, outer_radius, 10);
        for (unsigned int step = 0; step < 5; ++step)
        {
            for (auto &cell : triangulation.active_cell_iterators())
            {
                for (const auto v : cell->vertex_indices())
                {
                    const double distance_from_center =
                            center.distance(cell->vertex(v));
                    if (std::fabs(distance_from_center - inner_radius) <=
                        1e-6 * inner_radius)
                    {
                        cell->set_refine_flag();
                        break;
                    }
                }
            }
            triangulation.execute_coarsening_and_refinement();
        }
        std::ofstream out("grid-2.svg");
        GridOut grid_out;
        grid_out.write_svg(triangulation, out);
        std::cout << "Grid written to grid-2.svg" << std::endl;

}

*/
#endif //ION_PROPULSION_GRIDFORGE_H
