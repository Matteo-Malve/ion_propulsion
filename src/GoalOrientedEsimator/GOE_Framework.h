#ifndef GETPOT_GOE_FRAMEWORK_H
#define GETPOT_GOE_FRAMEWORK_H

#include "GOE_DualFunctional.h"

namespace GOE{
    using namespace dealii;
    template <int dim>


    struct Framework
    {
    public:
        using Evaluator     = Evaluation::EvaluationBase<dim>;
        using EvaluatorList = std::list<Evaluator *>;


        struct ProblemDescription
        {
            unsigned int primal_fe_degree;
            unsigned int dual_fe_degree;

            std::unique_ptr<const Data::SetUpBase<dim>> data;

            enum RefinementCriterion
            {
                dual_weighted_error_estimator,
                global_refinement,
                kelly_indicator,
                weighted_kelly_indicator
            };

            RefinementCriterion refinement_criterion;

            std::unique_ptr<const DualFunctional::DualFunctionalBase<dim>>
                    dual_functional;

            EvaluatorList evaluator_list;

            std::unique_ptr<const Function<dim>> kelly_weight;

            unsigned int max_degrees_of_freedom;

            ProblemDescription();
        };

        static void run(const ProblemDescription &descriptor);
    };


    template <int dim>
    Framework<dim>::ProblemDescription::ProblemDescription()
            : primal_fe_degree(1)
            , dual_fe_degree(2)
            , refinement_criterion(dual_weighted_error_estimator)
            , max_degrees_of_freedom(20000)
    {}



    template <int dim>
    void Framework<dim>::run(const ProblemDescription &descriptor)
    {
        Triangulation<dim> triangulation(
                Triangulation<dim>::smoothing_on_refinement);
        descriptor.data->create_coarse_grid(triangulation);

        const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
        const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
        const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
        const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);

        std::unique_ptr<LaplaceSolver::Base<dim>> solver;
        switch (descriptor.refinement_criterion)
        {
            case ProblemDescription::dual_weighted_error_estimator:
            {
                solver = std::make_unique<LaplaceSolver::WeightedResidual<dim>>(
                        triangulation,
                                primal_fe,
                                dual_fe,
                                quadrature,
                                face_quadrature,
                                descriptor.data->get_right_hand_side(),
                                descriptor.data->get_boundary_values(),
                                *descriptor.dual_functional);
                break;
            }

            case ProblemDescription::global_refinement:
            {
                solver = std::make_unique<LaplaceSolver::RefinementGlobal<dim>>(
                        triangulation,
                                primal_fe,
                                quadrature,
                                face_quadrature,
                                descriptor.data->get_right_hand_side(),
                                descriptor.data->get_boundary_values());
                break;
            }

            case ProblemDescription::kelly_indicator:
            {
                solver = std::make_unique<LaplaceSolver::RefinementKelly<dim>>(
                        triangulation,
                                primal_fe,
                                quadrature,
                                face_quadrature,
                                descriptor.data->get_right_hand_side(),
                                descriptor.data->get_boundary_values());
                break;
            }

            case ProblemDescription::weighted_kelly_indicator:
            {
                solver =
                        std::make_unique<LaplaceSolver::RefinementWeightedKelly<dim>>(
                                triangulation,
                                        primal_fe,
                                        quadrature,
                                        face_quadrature,
                                        descriptor.data->get_right_hand_side(),
                                        descriptor.data->get_boundary_values(),
                                        *descriptor.kelly_weight);
                break;
            }

            default:
                AssertThrow(false, ExcInternalError());
        }

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
}

#endif //GETPOT_GOE_FRAMEWORK_H
