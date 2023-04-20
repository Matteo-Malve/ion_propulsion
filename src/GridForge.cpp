#include "../include/GridForge.h"
//#include "GridForge.h"

// Grounbd mesh settings
// No 3,7,10,
// 1,2,4,5,12,13,14 unstructured
// 6,8 squadrata
// 9 Bordi
// 11 Vortici strani



template <int dim>
void CreateInitialGrid( Triangulation<dim> &mesh, const double mesh_height,
                 const double electrode_distance, const double wire_radius)
{
    Assert(dim == 2, ExcNotImplemented()); // Serve solo per dare errore se si prova con dim = 3
    /*
    const Point<2> top_left(0, mesh_height); // Coordinate angolo in alto a sinistra
    const Point<2> bottom_right(electrode_distance+wire_radius, 0.); // Coordinate angolo in basso a destra

    unsigned int a = (int) std::ceil( abs(top_left(0) - bottom_right(0))/wire_radius ); // Numero di suddivisioni lungo x
    unsigned int b = (int) std::ceil( abs(top_left(1) - bottom_right(1))/wire_radius ); // Numero di suddivisioni lungo y

    // Griglia rettangolare:
    GridGenerator::subdivided_hyper_rectangle( mesh, {a,b}, top_left, bottom_right);
    */

    // Check if ground mesh file is present and retrieve it
    struct stat sb;
    int is_present = stat("../gmsh_grids/custom_ground_mesh.msh",&sb);
    if(is_present==-1) {
        std::cerr << "Custom ground mesh NOT found!" << std::endl;
    }
    else if(is_present==0){
        std::cerr<<"Custom ground mesh found"<<std::endl;
        std::ifstream input_file("../gmsh_grids/custom_ground_mesh.msh");
        GridIn<dim>       grid_in;
        grid_in.attach_triangulation(mesh);
        grid_in.read_msh(input_file);
        std::cerr<<"Grid imported"<<std::endl;
    } else
        std::cerr<<"File not found nor not found. Anomaly.";

    // Initial refinement
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

    // Save initial mesh to file
    std::ofstream out("../mesh_storage/initial_mesh.vtu");
    GridOut       grid_out;
    //GridOutFlags::Vtu::serialize_triangulation = true;
    // void 	set_flags (const GridOutFlags::Vtu &flags)
    GridOutFlags::Vtu flags(true);
    grid_out.set_flags(flags);
    if(flags.serialize_triangulation==true)
        std::cout<<"GridOutFlags::Vtu::serialize_triangulation  IS  true"<<std::endl;
    else
        std::cout<<"GridOutFlags::Vtu::serialize_triangulation  IS  false"<<std::endl;
    grid_out.write_vtu(mesh, out);
    std::cout<<"Mesh written to vtu"<<std::endl;
}

#include <sys/stat.h>
template <int dim>
void CreateGrid( Triangulation<dim> &mesh, const double mesh_height,
                 const double electrode_distance, const double wire_radius){
    struct stat sb;
    int is_present = stat("../mesh_storage/initial_mesh.vtu",&sb);
    if(is_present==-1) {
        std::cout << "File NOT found: proceed to generate initial mesh" << std::endl;
        CreateInitialGrid<dim>(mesh, mesh_height, electrode_distance, wire_radius);
    }
    else if(is_present==0){
        std::cout<<"File found"<<std::endl;
        std::ifstream input_file("../mesh_storage/initial_mesh.vtu");
        GridIn<dim>       grid_in;
        grid_in.attach_triangulation(mesh);
        grid_in.read_vtu(input_file);
        std::cout<<"Grid imported"<<std::endl;
    } else
        std::cout<<"File not found nor not found. Anomaly.";
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