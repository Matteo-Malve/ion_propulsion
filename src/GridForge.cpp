#include "../include/GridForge.h"
//#include "GridForge.h"
template <int dim>
void CreateGrid( Triangulation<dim> &mesh, const double mesh_height,
                 const double electrode_distance, const double wire_radius)
{
    Assert(dim == 2, ExcNotImplemented()); // Serve solo per dare errore se si prova con dim = 3

    const Point<2> top_left(-electrode_distance-wire_radius, mesh_height); // Coordinate angolo in alto a sinistra
    const Point<2> bottom_right(electrode_distance+wire_radius, 0.); // Coordinate angolo in basso a destra
    //cout<<"Angolo alto a sinistra: ("<<top_left[0]<<","<<top_left[1]<<")"<<endl;
    //cout<<"Angolo basso a destra: ("<<bottom_right[0]<<","<<bottom_right[1]<<")"<<endl;

    unsigned int a = (int) std::ceil( (bottom_right(0) - top_left(0))/wire_radius ); // Numero di suddivisioni lungo x
    unsigned int b = (int) std::ceil( (top_left(1) - bottom_right(1))/wire_radius/2 ); // Numero di suddivisioni lungo y

    // Griglia rettangolare:
    GridGenerator::subdivided_hyper_rectangle( mesh, {a,b}, top_left, bottom_right);


    // Ciclo su tutte le facce per spostare quelle dove si trova l'elettrodo
    for (auto &face : mesh.active_face_iterators())
    {
        if (face->at_boundary()) // Se la faccia è sul bordo ...
        {
            const Point<2> c = face->center();
            if ( std::fabs(c[0]) <= wire_radius && c[1] < wire_radius*1e-3 ) // ... e il centro dista da 0,0 meno del raggio dell'emettitore
            {
                for (const auto i : face->vertex_indices()) // Ciclo sui vertici della cella
                {
                    Point<2> &v = face->vertex(i);
                    v(1) = std::max(0. , std::sqrt(pow(wire_radius,2) - pow(v(0),2) ) ); // impongo y = sqrt(r^2 - x^2)
                }
            }
        }
    }


}


template <int dim>
void SetManifoldsAndBoundaries(Triangulation<dim> &mesh, const double collector_height,
                               const double electrode_distance, const double wire_radius)
{
    Assert(dim == 2, ExcNotImplemented());

    // Boundaries
    const types::boundary_id collector_id = 2;
    const types::boundary_id emitter_id = 1;

    // Manifolds:
    const types::manifold_id     wire_id = 1;

    const Point<2>       center(0., 0.);
    SphericalManifold<dim> circle(center);
    mesh.set_manifold(wire_id, circle);

    for (auto &face : mesh.active_face_iterators())
    {
        if (face->at_boundary())
        {
            const Point<dim> c = face->center();

            if ( (c[1] > 0) && (c[1] <= collector_height) && (c[0] > electrode_distance) )
                face->set_boundary_id(collector_id);
            else if ( (c[1] > 0) && (std::fabs(c[0]) < 2*wire_radius) && (c[1] <= wire_radius*2)) {
                face->set_boundary_id(emitter_id);
                face->set_manifold_id(wire_id);
                //std::cout << " Set wire manifold in: " << c << std::endl;
            }
        }
    }
}

// #######################################
// Template initialization
// #######################################
template void CreateGrid( Triangulation<2> &mesh, const double mesh_height,
                          const double electrode_distance, const double wire_radius);
template void SetManifoldsAndBoundaries(Triangulation<2> &mesh, const double collector_height,
                                    const double electrode_distance, const double wire_radius);