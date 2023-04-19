#ifndef MESH_CC_PROBLEM_H
#define MESH_CC_PROBLEM_H

#include "includes&parameters_setup.h"
#include "HelperFunctions.h"
#include "GridForge.h"
#include <GetPot>
static GetPot datafile("../data_setup");

template <int dim>
class Problem
{
public:
    Problem() : fe(1), // elemeenti lineari (1) o quadratici (2)
                dof_handler(triangulation) {};

    // Ciclo globale: definizione matrici e condizioni al contorno
    void run();

private:
    // Creazione della griglia computazionale
    // Tutorial griglie: step-49 su deal.ii
    void create_mesh(const double mesh_height, const double electrode_distance,
                     const double wire_radius, const double collector_height);

    // Definizione delle matrici e dei vettori per la soluzione
    // Chiamata dopo ogni raffinimento della griglia per adeguare
    // la dimensione
    void setup_system();

    // Solver
    void solve();

    void refine_grid(const unsigned int min_grid_level, // minimo livello di raffinamento griglia
                     const unsigned int max_grid_level); // massimo livello di raffinemento griglia
    void output_results(const double wire_radius);

    Triangulation<dim> triangulation;   // step-1

    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;

    AffineConstraints<double> constraints;

    SparseMatrix<double> system_matrix;
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> laplace_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    Timer timer;

    unsigned int cycle = 0;
    const unsigned int Nmax = datafile("Numerics/FEM_cycles/Nmax",10); // massimo numero di cicli

    Vector<float> values;
    const float conv_tol = datafile("Numerics/FEM_cycles/global_tolerane",1e-4); // tolleranza globale
};


#endif //MESH_CC_PROBLEM_H
