#ifndef MESH_CC_GRIDFORGE_H
#define MESH_CC_GRIDFORGE_H

#include "includes&parameters_setup.h"



//  Funzione per creare la griglia
template <int dim>
void CreateGrid(Triangulation<dim> &mesh);
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


#endif //MESH_CC_GRIDFORGE_H
