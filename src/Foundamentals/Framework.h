#ifndef GETPOT_FRAMEWORK_H
#define GETPOT_FRAMEWORK_H

#include "../includes&parameters_setup.h"
#include "ErrorController.h"
#include "../Mesh/GridForge.h"

template <int dim>
struct Framework
{
public:
    struct ProblemDescription {
        unsigned int primal_fe_degree;
        unsigned int dual_fe_degree;
        std::unique_ptr<const DualFunctional::DualFunctionalBase <dim>> dual_functional;
        unsigned int max_degrees_of_freedom;
        Triangulation<dim> triangulation;
    }

    static void run(const ProblemDescription &descriptor);
};

template <int dim>
Framework<dim>::ProblemDescription::ProblemDescription()
        : primal_fe_degree(1)
        , dual_fe_degree(2)
        , max_degrees_of_freedom(20000)
{}

template <int dim>
void Framework<dim>::run(const ProblemDescription &descriptor)
{
    // MESH
    Triangulation<dim> triangulation;
    CreateGrid<dim>(triangulation);
    // FE
    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
    // QUADRATURE
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);
    // RHS
    Functions::Function<dim> rhs_function;

    //std::unique_ptr<Base<dim>> solver;


    solver = std::make_unique<ErrorController<dim>>(
            triangulation,
            primal_fe,
            dual_fe,
            quadrature,
            face_quadrature,
            rhs_function,
            descriptor.data->get_boundary_values(),
            *descriptor.dual_functional);


    for (unsigned int step = 0; true; ++step)
    {
        std::cout << "Refinement cycle: " << step << std::endl;

        solver->set_refinement_cycle(step);
        solver->solve_problem();
        solver->output_solution();

        std::cout << "   Number of degrees of freedom=" << solver->n_dofs()
                  << std::endl;

        for (const auto &evaluator : descriptor.evaluator_list)
        {
            evaluator->set_refinement_cycle(step);
            solver->postprocess(*evaluator);
        }


        if (solver->n_dofs() < descriptor.max_degrees_of_freedom)
            solver->refine_grid();
        else
            break;
    }

    std::cout << std::endl;
}


#endif //GETPOT_FRAMEWORK_H
