#include "globals.h"

namespace Evaluation {

template <int dim>
class EvaluationBase
{
public:
  virtual ~EvaluationBase() = default;
  void set_refinement_cycle(const unsigned int refinement_cycle);
  virtual double operator()(const DoFHandler<dim> &dof_handler,
                          const Vector<double>  &solution) const = 0;
protected:
  unsigned int refinement_cycle;
};

template <int dim>
void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step)
{
  refinement_cycle = step;
}

template <int dim>
class PointValueEvaluation : public EvaluationBase<dim>
{
public:
  PointValueEvaluation(const Point<dim> &evaluation_point);
  virtual double operator()(const DoFHandler<dim> &dof_handler,
                          const Vector<double>  &solution) const override;
  DeclException1(
    ExcEvaluationPointNotFound,
    Point<dim>,
    << "The evaluation point " << arg1
    << " was not found among the vertices of the present grid.");
private:
  const Point<dim> evaluation_point;
};

template <int dim>
PointValueEvaluation<dim>::PointValueEvaluation(
  const Point<dim> &evaluation_point)
  : evaluation_point(evaluation_point)
{}

template <int dim>
double PointValueEvaluation<dim>::operator()(const DoFHandler<dim> &dof_handler, 
                                          const Vector<double>  &solution) const
{
  // Iterate over all active cells to find the vertex
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
    {
      if (cell->vertex(v).distance(evaluation_point) < 1e-12) // Match vertex point
      {
        // Get the global DoF index for the vertex
        unsigned int vertex_dof_index = cell->vertex_dof_index(v, 0);
        
        // Return the solution value at that vertex
        return solution(vertex_dof_index);
      }
    }
  }
  //throw std::runtime_error("Vertex not found in the mesh.");
  cout<<"      Vertex not found in the mesh."<<endl;
  return VectorTools::point_value(dof_handler, solution, evaluation_point);
}



template <int dim>
class PointYDerivativeEvaluation : public EvaluationBase<dim>
{
public:
  PointYDerivativeEvaluation(const Point<dim> &evaluation_point);
  virtual double operator()(const DoFHandler<dim> &dof_handler,
                          const Vector<double>  &solution) const override;
  DeclException1(
    ExcEvaluationPointNotFound,
    Point<dim>,
    << "The evaluation point " << arg1
    << " was not found among the vertices of the present grid.");
private:
  const Point<dim> evaluation_point;
};

template <int dim>
PointYDerivativeEvaluation<dim>::PointYDerivativeEvaluation(
  const Point<dim> &evaluation_point)
  : evaluation_point(evaluation_point)
{}

/*template <int dim>
double PointYDerivativeEvaluation<dim>::operator()(
  const DoFHandler<dim> &dof_handler,
  const Vector<double>  &solution) const
{
  double point_derivative = 0;
  const QTrapezoid<dim>       vertex_quadrature;
  FEValues<dim>               fe_values(dof_handler.get_fe(),
                          vertex_quadrature,
                          update_gradients | update_quadrature_points);
  std::vector<Tensor<1, dim>> solution_gradients(vertex_quadrature.size());
  unsigned int evaluation_point_hits = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    for (const auto vertex : cell->vertex_indices())
      if (cell->vertex(vertex) == evaluation_point)
        {
          fe_values.reinit(cell);
          fe_values.get_function_gradients(solution, solution_gradients);
          unsigned int q_point = 0;
          for (; q_point < solution_gradients.size(); ++q_point)
            if (fe_values.quadrature_point(q_point) == evaluation_point)
              break;
          Assert(q_point < solution_gradients.size(), ExcInternalError());
          point_derivative += solution_gradients[q_point][0];
          ++evaluation_point_hits;
break; }
  AssertThrow(evaluation_point_hits > 0,
              ExcEvaluationPointNotFound(evaluation_point));
  point_derivative /= evaluation_point_hits;
  std::cout << "   Point x-derivative=" << point_derivative << std::endl;
}*/


template <int dim>
double PointYDerivativeEvaluation<dim>::operator()(
  const DoFHandler<dim> &dof_handler,
  const Vector<double>  &solution) const
{
  double point_derivative = 0;
  const QTrapezoid<dim>       vertex_quadrature;
  FEValues<dim>               fe_values(dof_handler.get_fe(),
                          vertex_quadrature,
                          update_gradients | update_quadrature_points);
  std::vector<Tensor<1, dim>> solution_gradients(vertex_quadrature.size());

  // Variables to track the closest point if exact match is not found
  double min_distance = std::numeric_limits<double>::max();
  typename DoFHandler<dim>::active_cell_iterator nearest_cell;
  unsigned int nearest_vertex = 0;
  bool exact_match_found = false;

  // Placeholder for the final point to be used in evaluation
  Point<dim> adjusted_evaluation_point = evaluation_point;

  // First loop to find the exact evaluation point or the closest one if it doesn't exist
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    for (const auto vertex : cell->vertex_indices()) {
      double distance = cell->vertex(vertex).distance(evaluation_point);

      if (distance < 1e-6 * cell->diameter()) {
        // Exact point found
        exact_match_found = true;
        break;
      }

      // Track the nearest vertex
      if (distance < min_distance) {
        min_distance = distance;
        nearest_cell = cell;
        nearest_vertex = vertex;
      }
    }
    if (exact_match_found) break;
  }

  // Set the evaluation point to the nearest vertex if no exact match was found
  if (!exact_match_found) {
    adjusted_evaluation_point = nearest_cell->vertex(nearest_vertex);
    cout<<"         Adjusted the evaluation point"<<endl;

  }


  unsigned int evaluation_point_hits = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    for (const auto vertex : cell->vertex_indices())
      if (cell->vertex(vertex).distance(adjusted_evaluation_point) < 1e-6 * cell->diameter())
        {
          fe_values.reinit(cell);
          fe_values.get_function_gradients(solution, solution_gradients);
          unsigned int q_point = 0;
          for (; q_point < solution_gradients.size(); ++q_point)
            if (fe_values.quadrature_point(q_point).distance(adjusted_evaluation_point) < 1e-6 * cell->diameter())
              break;

          Assert(q_point < solution_gradients.size(), ExcInternalError());
          point_derivative += solution_gradients[q_point][1];
          ++evaluation_point_hits;
          
          break; 
        }
  AssertThrow(evaluation_point_hits > 0,
              ExcEvaluationPointNotFound(adjusted_evaluation_point));
  point_derivative /= evaluation_point_hits;
  return point_derivative;
}


}