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


#endif //ION_PROPULSION_GRIDFORGE_H
