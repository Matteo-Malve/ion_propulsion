  VectorTools::project_boundary_values(MappingQ1<dim>(), primal_dof_handler, {{1, &exact_solution_function}}, QGauss<dim-1>(4), emitter_boundary_values);



template <int dim>
double get_solution_at_vertex(const DoFHandler<dim> &dof_handler,
                              const Vector<double> &solution,
                              const Point<dim> &vertex_point)
{
  // Iterate over all active cells to find the vertex
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
    {
      if (cell->vertex(v).distance(vertex_point) < 1e-12) // Match vertex point
      {
        // Get the global DoF index for the vertex
        unsigned int vertex_dof_index = cell->vertex_dof_index(v, 0);
        
        // Return the solution value at that vertex
        return solution(vertex_dof_index);
      }
    }
  }
  //throw std::runtime_error("Vertex not found in the mesh.");
  cout<<"Vertex not found in the mesh."<<endl;
  return VectorTools::point_value(dof_handler, solution, vertex_point);
}