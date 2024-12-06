#include "Data.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace Data{

    // ------------------------------------------------------
    // SetUp
    // ------------------------------------------------------

    template <class Traits, int dim>
    const Function<dim> &SetUp<Traits, dim>::get_boundary_values() const
    {
      return boundary_values;
    }

    template <class Traits, int dim>
    const Function<dim> &SetUp<Traits, dim>::get_right_hand_side() const
    {
      return right_hand_side;
    }


    template <class Traits, int dim>
    void SetUp<Traits, dim>::create_coarse_grid(
      Triangulation<dim> &coarse_grid) const
    {
      Traits::create_coarse_grid(coarse_grid);
    }

    // ------------------------------------------------------
    // CurvedRidges
    // ------------------------------------------------------

    template <int dim>
    double CurvedRidges<dim>::BoundaryValues::value(
      const Point<dim> &p,
      const unsigned int /*component*/) const
    {
      double q = p(0);
      for (unsigned int i = 1; i < dim; ++i)
        q += std::sin(10 * p(i) + 5 * p(0) * p(0));
      const double exponential = std::exp(q);
      return exponential;
    }

    template <int dim>
    double CurvedRidges<dim>::RightHandSide::value(
      const Point<dim> &p,
      const unsigned int /*component*/) const
    {
      double q = p(0);
      for (unsigned int i = 1; i < dim; ++i)
        q += std::sin(10 * p(i) + 5 * p(0) * p(0));
      const double u  = std::exp(q);
      double       t1 = 1, t2 = 0, t3 = 0;
      for (unsigned int i = 1; i < dim; ++i)
      {
        t1 += std::cos(10 * p(i) + 5 * p(0) * p(0)) * 10 * p(0);
        t2 += 10 * std::cos(10 * p(i) + 5 * p(0) * p(0)) -
              100 * std::sin(10 * p(i) + 5 * p(0) * p(0)) * p(0) * p(0);
        t3 += 100 * std::cos(10 * p(i) + 5 * p(0) * p(0)) *
                std::cos(10 * p(i) + 5 * p(0) * p(0)) -
              100 * std::sin(10 * p(i) + 5 * p(0) * p(0));
      }
      t1 = t1 * t1;

      return -u * (t1 + t2 + t3);
    }

    template <int dim>
    void CurvedRidges<dim>::create_coarse_grid(Triangulation<dim> &coarse_grid)
    {
      GridGenerator::hyper_cube(coarse_grid, -1, 1);
      coarse_grid.refine_global(2);
    }

    // ------------------------------------------------------
    // Exercise_2_3
    // ------------------------------------------------------

    template <>
    void Exercise_2_3<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      const unsigned int dim = 2;

      const std::vector<Point<2>> vertices = {
        {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0}, //
        {-1.0, -0.5}, {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5}, //
        {-1.0, +0.0}, {-0.5, +0.0}, {+0.5, +0.0}, {+1.0, +0.0},               //
        {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5}, {+1.0, +0.5}, //
        {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}};

      // Next, we have to define the cells and the vertices they contain.
      const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>>
        cell_vertices = {{{0, 1, 5, 6}},
                         {{1, 2, 6, 7}},
                         {{2, 3, 7, 8}},
                         {{3, 4, 8, 9}},
                         {{5, 6, 10, 11}},
                         {{8, 9, 12, 13}},
                         {{10, 11, 14, 15}},
                         {{12, 13, 17, 18}},
                         {{14, 15, 19, 20}},
                         {{15, 16, 20, 21}},
                         {{16, 17, 21, 22}},
                         {{17, 18, 22, 23}}};

      const unsigned int n_cells = cell_vertices.size();

      // Again, we generate a C++ vector type from this, but this time by
      // looping over the cells (yes, this is boring). Additionally, we set
      // the material indicator to zero for all the cells:
      std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
      for (unsigned int i = 0; i < n_cells; ++i)
        {
          for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
            cells[i].vertices[j] = cell_vertices[i][j];
          cells[i].material_id = 0;
        }

      // Finally pass all this information to the library to generate a
      // triangulation. The last parameter may be used to pass information
      // about non-zero boundary indicators at certain faces of the
      // triangulation to the library, but we don't want that here, so we give
      // an empty object:
      coarse_grid.create_triangulation(vertices, cells, SubCellData());

      // And since we want that the evaluation point (3/4,3/4) in this example
      // is a grid point, we refine once globally:
      coarse_grid.refine_global(1);
    }
  // Template instantiation
  template struct SetUpBase<2>;
  template struct SetUp<IonPropulsion::Data::Exercise_2_3<2>, 2>;
  template struct CurvedRidges<2>;
  template struct Exercise_2_3<2>;
  } // namespace Data
} // namespace IonPropulsion