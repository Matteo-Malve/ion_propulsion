#ifndef GETPOT_FRAMEWORK_H
#define GETPOT_FRAMEWORK_H

#include "../includes&parameters_setup.h"
#include "ErrorController.h"
#include "../Mesh/GridForge.h"
#include "Base.h"
#include "DualFunctional.h"

// Forward declaration of ErrorController class
template <int dim>
class ErrorController;


template<int dim>
struct ProblemDescription {
    unsigned int primal_fe_degree;
    unsigned int dual_fe_degree;
    std::unique_ptr<const DualFunctionalBase<dim>> dual_functional;
    unsigned int max_degrees_of_freedom;
    Triangulation<dim> triangulation;
    // Tripletta infame
    Function<dim> &      rhs_function;        // Not sure
    //const Function<dim> &      boundary_values;     // Not sure
    // Constructor
    ProblemDescription(Function<dim>&);
};

template<int dim>
ProblemDescription<dim>::ProblemDescription(Function<dim>& rhs_func)
        : primal_fe_degree(1)
        , dual_fe_degree(2)
        , max_degrees_of_freedom(50000),
          rhs_function(rhs_func) {}


template <int dim,unsigned int grid_option>
struct Framework{
public:
    static void run(const ProblemDescription<dim> &descriptor);
};

template <int dim, unsigned int grid_option>
void Framework<dim,grid_option>::run(const ProblemDescription<dim> &descriptor)

{
    // MESH
    Triangulation<dim> triangulation(Triangulation<dim>::smoothing_on_refinement);
    switch(grid_option){
        case 1:
            CreateGrid<dim>(triangulation);
        case 2:
            LoadSecondGrid<dim>(triangulation);
        default:
            cout<<"ERROR2: Invalid grid_option id";
    }



    // FE
    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
    // QUADRATURE
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);

    // SOLVER = pointer to BASE
    std::unique_ptr<Base<dim>> solver;
    solver = std::make_unique<ErrorController<dim>>(
            triangulation,
            primal_fe,
            dual_fe,
            quadrature,
            face_quadrature,
            descriptor.rhs_function,    // Not sure
            //descriptor.boundary_values, // Not sure
            *descriptor.dual_functional);


    // REFINEMENT LOOP
    for (unsigned int step = 0; true; ++step)
    {
        std::cout << "[Framework]Refinement cycle: " << step << std::endl;

        solver->set_refinement_cycle(step);
        solver->solve_problem();    // <--- Problema qui
        solver->output_solution();


        std::cout << "   [Framework]Number of degrees of freedom=" << solver->n_dofs()
                  << std::endl;

        if (solver->n_dofs() < descriptor.max_degrees_of_freedom) {
            cout << "   [Framework]Prepare call to refine grid" << endl;
            solver->refine_grid();
        }
        else
            break;
    }

    std::cout << std::endl;
}

/*
// Specialization for grid:option=1
template<int dim>
void Framework<dim,1>::run(const ProblemDescription<dim> &descriptor)
{
    // MESH
    Triangulation<dim> triangulation(Triangulation<dim>::smoothing_on_refinement);
    CreateGrid<dim>(triangulation);
    //LoadSecondGrid<dim>(triangulation);
    
    // FE
    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
    // QUADRATURE
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);

    // SOLVER = pointer to BASE
    std::unique_ptr<Base<dim>> solver;
    solver = std::make_unique<ErrorController<dim>>(
            triangulation,
            primal_fe,
            dual_fe,
            quadrature,
            face_quadrature,
            descriptor.rhs_function,    // Not sure
            //descriptor.boundary_values, // Not sure
            *descriptor.dual_functional);


    // REFINEMENT LOOP
    for (unsigned int step = 0; true; ++step)
    {
        std::cout << "[Framework]Refinement cycle: " << step << std::endl;

        solver->set_refinement_cycle(step);
        solver->solve_problem();    // <--- Problema qui
        solver->output_solution();


        std::cout << "   [Framework]Number of degrees of freedom=" << solver->n_dofs()
                  << std::endl;

        if (solver->n_dofs() < descriptor.max_degrees_of_freedom) {
            cout << "   [Framework]Prepare call to refine grid" << endl;
            solver->refine_grid();
        }
        else
            break;
    }

    std::cout << std::endl;
}

// Specialization for grid_option=2
template<int dim>
void Framework<dim,2>::run(const ProblemDescription<dim> &descriptor)
{
    // MESH
    Triangulation<dim> triangulation(Triangulation<dim>::smoothing_on_refinement);
    //CreateGrid<dim>(triangulation);
    LoadSecondGrid<dim>(triangulation);

    // FE
    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
    // QUADRATURE
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);

    // SOLVER = pointer to BASE
    std::unique_ptr<Base<dim>> solver;
    solver = std::make_unique<ErrorController<dim>>(
            triangulation,
            primal_fe,
            dual_fe,
            quadrature,
            face_quadrature,
            descriptor.rhs_function,    // Not sure
            //descriptor.boundary_values, // Not sure
            *descriptor.dual_functional);


    // REFINEMENT LOOP
    for (unsigned int step = 0; true; ++step)
    {
        std::cout << "[Framework]Refinement cycle: " << step << std::endl;

        solver->set_refinement_cycle(step);
        solver->solve_problem();    // <--- Problema qui
        solver->output_solution();


        std::cout << "   [Framework]Number of degrees of freedom=" << solver->n_dofs()
                  << std::endl;

        if (solver->n_dofs() < descriptor.max_degrees_of_freedom) {
            cout << "   [Framework]Prepare call to refine grid" << endl;
            solver->refine_grid();
        }
        else
            break;
    }

    std::cout << std::endl;
}
*/

#endif //GETPOT_FRAMEWORK_H
