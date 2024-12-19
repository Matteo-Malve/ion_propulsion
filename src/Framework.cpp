#include "Framework.h"

namespace IonPropulsion{
  using namespace dealii;

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
    // First create a triangulation from the given data object,
    Triangulation<dim> triangulation(
      Triangulation<dim>::smoothing_on_refinement);
    descriptor.data->create_coarse_grid(triangulation);

    // then a set of finite elements and appropriate quadrature formula:
    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);

    // Next, select one of the classes implementing different refinement
    // criteria.
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
      //solver->compute_flux();
      solver->update_convergence_table();
      solver->output_solution();

      std::cout << "   Number of degrees of freedom=" << solver->n_dofs()
                << std::endl;

      for (const auto &evaluator : descriptor.evaluator_list)
        {
          evaluator->set_refinement_cycle(step);
          solver->postprocess(*evaluator);
        }

      unsigned int DoFs_before_refinement = solver->n_dofs();
      solver->refine_grid();
      solver->print_convergence_table();
      CSVLogger::getInstance().flushRow();

      if (DoFs_before_refinement > descriptor.max_degrees_of_freedom)
        break;

    }

  // Clean up the screen after the loop has run:
  std::cout << std::endl;
}

  // Template instantiation
  template struct Framework<2>;
} // namespace IonPropulsion