#include <iostream>
#include "../src/includes&parameters_setup.h"


int main()
{
    // BUILD TRIANGULATION
    Triangulation<2> triangulation;
    const Point<2> center(0, 0);
    const double inner_radius = 0.025, outer_radius = 1.025;
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


    // SET BOUNDARY
    const types::boundary_id collector_id = 2;
    const types::boundary_id emitter_id = 1;

    for (auto &face : triangulation.active_face_iterators())
    {
        if (face->at_boundary())
        {
            const double distance_from_center = center.distance(face->center());
            if (distance_from_center <= 1.02 * inner_radius)
                face->set_boundary_id(emitter_id);
            else if (distance_from_center >= 0.98*outer_radius) {
                face->set_boundary_id(collector_id);
            }
        }
    }

    // OUTPUT
    std::ofstream out("grid-2.svg");
    GridOut grid_out;
    grid_out.write_svg(triangulation, out);
    std::cout << "Grid written to grid-2.svg" << std::endl;

    std::ofstream out2("circular_mesh.vtu");
    GridOut       grid_out2;
    GridOutFlags::Vtu flags(true);
    grid_out2.set_flags(flags);
    cout<<endl<<"Saving constructed mesh to file:"<<endl;
    if(flags.serialize_triangulation==true)
        std::cout<<" GridOutFlags::Vtu::serialize_triangulation  IS  true"<<std::endl;
    else
        std::cout<<" GridOutFlags::Vtu::serialize_triangulation  IS  false"<<std::endl;
    grid_out2.write_vtu(triangulation, out2);
    std::cout<<" Mesh written to vtu"<<endl<<endl;
    return 0;
}