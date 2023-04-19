#include "../include/GridForge.h"
//#include "GridForge.h"

template <int dim>
void CreateInitialGrid( Triangulation<dim> &mesh, const double mesh_height,
                 const double electrode_distance, const double wire_radius)
{
    Assert(dim == 2, ExcNotImplemented()); // Serve solo per dare errore se si prova con dim = 3

    const Point<2> top_left(0, mesh_height); // Coordinate angolo in alto a sinistra
    const Point<2> bottom_right(electrode_distance+wire_radius, 0.); // Coordinate angolo in basso a destra

    unsigned int a = (int) std::ceil( abs(top_left(0) - bottom_right(0))/wire_radius ); // Numero di suddivisioni lungo x
    unsigned int b = (int) std::ceil( abs(top_left(1) - bottom_right(1))/wire_radius ); // Numero di suddivisioni lungo y

    // Griglia rettangolare:
    GridGenerator::subdivided_hyper_rectangle( mesh, {a,b}, top_left, bottom_right);

    const Point<2> center(0, 0);
    for (unsigned int step = 0; step < 5; ++step)
    {
        for (auto &cell : mesh.active_cell_iterators())
        {
            const double distance_from_center = center.distance(cell->center());
            if(distance_from_center<0.5)
                cell->set_refine_flag();
            else if(step<4 && distance_from_center<0.75)
                cell->set_refine_flag();
            else if(step<3 && distance_from_center<1)
                cell->set_refine_flag();
            else if(step<2 && distance_from_center<1.25)
                cell->set_refine_flag();
            }
        mesh.execute_coarsening_and_refinement();
    }
    std::ofstream out("../mesh_storage/initial_mesh.vtu");
    GridOut       grid_out;
    grid_out.write_vtu(mesh, out);
    std::cout<<"Mesh written to vtu"<<std::endl;
}

template <int dim>
void CreateGrid( Triangulation<dim> &mesh, const double mesh_height,
                 const double electrode_distance, const double wire_radius){
    CreateInitialGrid<dim>(mesh,mesh_height,electrode_distance,wire_radius);
    std::cout<<"Tutto bene"<<std::endl;
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