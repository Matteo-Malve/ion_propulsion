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
    unsigned int max_number_of_refinement_cycles;
    Triangulation<dim> triangulation;
    Function<dim> &      rhs_function;
    // Constructor
    ProblemDescription(Function<dim>&);
};

template<int dim>
ProblemDescription<dim>::ProblemDescription(Function<dim>& rhs_func)
        : primal_fe_degree(1)
        , dual_fe_degree(2)
        , max_degrees_of_freedom(20000)
        , max_number_of_refinement_cycles(0),
          rhs_function(rhs_func) {}

template<int dim>
void framework_run(const ProblemDescription<dim> &descriptor,unsigned int grid_option){
        Triangulation<dim> triangulation(Triangulation<dim>::smoothing_on_refinement);
    // IMPORT correct grid
        if(grid_option<1.5)
            CreateGrid<dim>(triangulation);
        else
            LoadSecondGrid<dim>(triangulation);
    // FE
    const FE_Q <dim> primal_fe(descriptor.primal_fe_degree);
    const FE_Q <dim> dual_fe(descriptor.dual_fe_degree);
    // QUADRATURE
    const QGauss <dim> quadrature(descriptor.dual_fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);

    // SOLVER = pointer to BASE
    std::unique_ptr <Base<dim>> solver;
    solver = std::make_unique<ErrorController<dim>>(
            triangulation,
            primal_fe,
            dual_fe,
            quadrature,
            face_quadrature,
            descriptor.rhs_function,
            *descriptor.dual_functional);


    // REFINEMENT LOOP: ruled both by -Max Nb. of iterations- and by a -Max Nb. of DoF-
    unsigned int step = 0;
    while(step<descriptor.max_number_of_refinement_cycles && solver->n_dofs() < descriptor.max_degrees_of_freedom){
        std::cout << "[Framework]Refinement cycle: " << step << std::endl;

        // Core calls
        solver->set_refinement_cycle(step);
        solver->solve_problem();
        solver->output_solution();
        std::cout << "   [Framework]Number of degrees of freedom=" << solver->n_dofs() << std::endl;

        // Refine grid
        cout << "   [Framework]Prepare call to refine grid" << endl;
        solver->refine_grid(step);

        // Update step
        step++;
    }

    std::cout << std::endl;
}


#endif //GETPOT_FRAMEWORK_H
