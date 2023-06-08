#ifndef MESH_CC_PROBLEM_H
#define MESH_CC_PROBLEM_H

#include "../includes&parameters_setup.h"
#include "HelperFunctions.h"
#include "../Mesh/GridForge.h"
#include "DirichletBoundaryValuesDX.h"
#include "../Evaluation.h"
static GetPot datafile("../data_setup");

template<int dim>
class Solverbase{
public:
    Solverbase(const unsigned int fe_order): fe(fe_order), dof_handler(triangulation) {};
    virtual void run();
    
    Triangulation<dim> triangulation;
    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;
    Vector<double> solution;
    

protected:
    void create_mesh();
    void setup_system();
    void apply_boundary_conditions();
    virtual void assemble_rhs(){
        VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), system_rhs);
    }
    virtual void save_finest_mesh(){};

    // Solver
    void solve();

    void refine_grid(const unsigned int min_grid_level, // minimo livello di raffinamento griglia
                     const unsigned int max_grid_level); // massimo livello di raffinemento griglia
    void output_results();

    
    Triangulation<dim> finest_mesh;

    

    AffineConstraints<double> constraints;

    SparseMatrix<double> system_matrix;
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> laplace_matrix;

    
    Vector<double> system_rhs;

    Timer timer;

    unsigned int cycle = 0;
    unsigned int Nmax = datafile("Numerics/FEM_cycles/Nmax",10);

    Vector<float> values;
    const float conv_tol = datafile("Numerics/FEM_cycles/global_tolerane",1e-4); // tolleranza globale

    double wire_radius = datafile("Mesh/wire_radius",0.025);

};


#endif //MESH_CC_PROBLEM_H
