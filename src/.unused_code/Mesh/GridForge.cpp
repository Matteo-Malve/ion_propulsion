#include "GridForge.h"
#include <sys/stat.h>

template <int dim>
void CreateInitialGrid(Triangulation<dim> &mesh)
{
    Assert(dim == 2, ExcNotImplemented());

    // Check if ground mesh file is present and retrieve it
    cout<<"[GridForge::CreateInitialGrid]Looking for Custom Ground Mesh built with Gmsh"<<endl;
    struct stat sb;
    int is_present = stat("../gmsh_grids/input_mesh.msh",&sb);
    if(is_present==-1) {
        std::cerr << "   Mesh NOT found!" << std::endl;
    }
    else if(is_present==0){
        std::cerr<<"   Mesh found"<<std::endl<<"   Prepare import"<<endl;
        std::ifstream input_file("../gmsh_grids/input_mesh.msh");
        GridIn<dim>       grid_in;
        grid_in.attach_triangulation(mesh);
        grid_in.read_msh(input_file);
        std::cerr<<"   Grid imported"<<std::endl;
    } else
        std::cerr<<"   File not found nor not found. Anomaly."<<endl;

    // Save initial mesh to file
    std::ofstream out("../mesh_storage/initial_mesh.vtu");
    GridOut       grid_out;
    GridOutFlags::Vtu flags(true);
    grid_out.set_flags(flags);
    cout<<endl<<"Saving constructed mesh to file:"<<endl;
    if(flags.serialize_triangulation==true)
        std::cout<<"   GridOutFlags::Vtu::serialize_triangulation  IS  true"<<std::endl;
    else
        std::cout<<"   GridOutFlags::Vtu::serialize_triangulation  IS  false"<<std::endl;
    grid_out.write_vtu(mesh, out);
    std::cout<<"   Mesh written to .vtu"<<endl<<endl;

    // NOTE:  To clean the mesh_storage folder and import new grids, remember to type the command:
    //        > make cleanmesh
    //        before the usual
    //        > make run
}


// If a grid is already stored, it won't be necessary to import a new one
template <int dim>
void CreateGrid(Triangulation<dim> &mesh){
    struct stat sb;
    int is_present = stat("../mesh_storage/initial_mesh.vtu",&sb);
    cout<<endl<<"[GridForge::CreateGrid]Looking for an already existent mesh:"<<endl;
    if(is_present==-1) {
        std::cout << "   [GridForge::CreateGrid]File NOT found: proceed to generate initial mesh" << std::endl;
        CreateInitialGrid<dim>(mesh);
    }
    else if(is_present==0){
        std::cout<<"   [GridForge::CreateGrid]File found"<<std::endl;
        std::cout<<"   [GridForge::CreateGrid]Prepare import"<<std::endl;
        std::ifstream input_file("../mesh_storage/initial_mesh.vtu");
        GridIn<dim>       grid_in;
        grid_in.attach_triangulation(mesh);
        grid_in.read_vtu(input_file);
        std::cout<<"   [GridForge::CreateGrid]Grid imported"<<std::endl;
    } else
        std::cout<<"   [GridForge::CreateGrid]File not found nor not found. Anomaly."<<endl;
}


template <int dim>
void LoadSecondGrid(Triangulation<dim> &mesh){
    struct stat sb;
    int is_present = stat("../gmsh_grids/circular_mesh.msh",&sb);
    cout<<endl<<"[GridForge::LoadSecondGrid]Looking for an already existent mesh:"<<endl;
    if(is_present==-1) {
        std::cout << "   [GridForge::LoadSecondGrid]File NOT found: proceed to generate initial mesh" << std::endl;
    }
    else if(is_present==0){
        std::cout<<"   [GridForge::LoadSecondGrid]File found"<<std::endl;
        std::cout<<"   [GridForge::LoadSecondGrid]Prepare import"<<std::endl;
        std::ifstream input_file("../gmsh_grids/circular_mesh.msh");
        GridIn<dim>       grid_in;
        grid_in.attach_triangulation(mesh);
        grid_in.read_msh(input_file);
        std::cout<<"   [GridForge::LoadSecondGrid]Grid imported"<<std::endl;

    } else
        std::cout<<"   [GridForge::LoadSecondGrid]File not found nor not found. Anomaly."<<endl;
}



// #######################################
// Template initialization
// #######################################
template void CreateGrid( Triangulation<2> &mesh);
template void LoadSecondGrid( Triangulation<2> &mesh);