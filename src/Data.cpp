#include "Data.h"
#include <deal.II/grid/grid_in.h>

namespace IonPropulsion {
  using namespace dealii;
  using std::endl;
  namespace Data {
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
      static const typename Traits::RightHandSide rhs; // Lazy initialization
      return right_hand_side;
    }

    template <class Traits, int dim>
    const Function<dim> &SetUp<Traits, dim>::get_exact_solution() const
    {
      return exact_solution;
    }


    template <class Traits, int dim>
    void SetUp<Traits, dim>::create_coarse_grid(
      Triangulation<dim> &coarse_grid) const
    {
      Traits::create_coarse_grid(coarse_grid);
    }

    // ------------------------------------------------------
    // SetupNone
    // ------------------------------------------------------

    template <int dim,int sign>
    class CollectorGeometry : public ChartManifold<dim, dim, dim-1>
    {
    public:
      virtual Point<dim-1> pull_back(const Point<dim> &space_point) const override;

      virtual Point<dim> push_forward(const Point<dim-1> &chart_point) const override;

      virtual std::unique_ptr<Manifold<dim, dim>> clone() const override;

    };

    // HELPER FUNCTION
    double get_collector_height(const double &p);
    Tensor<1,2> get_emitter_normal(const Point<2> a);

    template <int dim,int sign>
    std::unique_ptr<Manifold<dim, dim>> CollectorGeometry<dim,sign>::clone() const
    {
      return std::make_unique<CollectorGeometry<dim,sign>>();
    }

    //--------------------------------------------------------------------------------------------------------------------------------------------

    template <int dim,int sign>
    Point<dim> CollectorGeometry<dim,sign>::push_forward(const Point<dim-1>  &x) const
    {
      const double y = sign*get_collector_height( x[0] );

      Point<dim> p;
      p[0] = x[0]; p[1] = y;

      if (dim == 3) {
        p[2] = x[1];
      }

      return p;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    template <int dim, int sign>
    Point<dim-1>  CollectorGeometry<dim, sign>::pull_back(const Point<dim> &p) const
    {
      Point<dim-1> x;
      x[0] = p[0];

      if (dim == 3) {
        x[1] = p[2];
      }

      return x;
    }

    // ############################ HELPER FUNCTIONS ##################################################################################

    double X = -2.5e-5;
    double g = 0.02;
    double collector_length = 0.10;

    double get_collector_height(const double &X)
    {
      const double x = (X-g)/collector_length;
      double y = 0;

      if ( abs(x-1.) > 1e-12 && abs(x) > 1e-12 ) {
        double a0 = 0.2969;
        double a1 = -0.126;
        double a2 = -0.3516;
        double a3 = 0.2843;
        double a4 = -0.1036; // or -0.1015 for an open trailing edge
        double t = 0.5; // Last 2 digits of the NACA divided by 20

        y = t*( a0 * sqrt(x) + a1 * x + a2 * pow(x,2.0) + a3 * pow(x,3.0) + a4 * pow(x,4.0) );
      }

      return y * collector_length;
    }

    //-------------------------------------------------------------------------------------------------------------------------------------------------

    Tensor<1,2> get_emitter_normal(const Point<2> a) {

      Tensor<1,2> normal;

      normal[0] = a[0] - X;
      normal[1] = a[1];

      const double norm = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);

      return normal/norm;
    }

    //-------------------------------------------------------------------------------------------------------------------------------------------------


    template <int dim>
    void SetupNone<dim>::create_coarse_grid(Triangulation<dim> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      const types::manifold_id emitter = 1;
      const Point<2> center(X,0.0);
      SphericalManifold<2> emitter_manifold(center);

      const types::manifold_id upper_collector = 2;
      CollectorGeometry<2,+1> upper_collector_manifold;

      const types::manifold_id bottom_collector = 3;
      CollectorGeometry<2,-1> bottom_collector_manifold;

      coarse_grid.set_all_manifold_ids_on_boundary(1, emitter);
      coarse_grid.set_manifold(emitter, emitter_manifold);
      coarse_grid.set_all_manifold_ids_on_boundary(2, upper_collector);
      coarse_grid.set_manifold(upper_collector, upper_collector_manifold);
      coarse_grid.set_all_manifold_ids_on_boundary(3, bottom_collector);
      coarse_grid.set_manifold(bottom_collector, bottom_collector_manifold);

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

      return - eps_0 * eps_r * u * (t1 + t2 + t3);  //  TODO not sure about eps
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

      for(auto cell : coarse_grid.active_cell_iterators())
        for (unsigned int face = 0; face < GeometryInfo<2>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(9);
      // And since we want that the evaluation point (3/4,3/4) in this example
      // is a grid point, we refine once globally:
      coarse_grid.refine_global(1);
    }

    // ------------------------------------------------------
    // Rectangle_1_99
    // ------------------------------------------------------

    template <>
    void Rectangle_1_99<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      //const std::string path_to_mesh = "../mesh/TinyStep14_1_99.msh";
      //const std::string path_to_mesh = "../mesh/TinyStep14_deFalco.msh";
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      double l = 0.0001;

      const Point<2> center(0.0, 0.0);
      double D = 4*l;
      for (unsigned int i = 0; i < NUM_CONCENTRIC_REF; ++i) {
        Vector<float> criteria(coarse_grid.n_active_cells());
        unsigned int ctr = 0;

        D *= 0.999;

        for (auto &cell : coarse_grid.active_cell_iterators()) {
          Point<2> closest_vertex_p = cell->center();
          double min_distance = closest_vertex_p.distance(center);
          for (const auto vertex : cell->vertex_indices()) {
            auto & vertex_p = cell->vertex(vertex);
            if (vertex_p.distance(center)<min_distance) {
              min_distance = vertex_p.distance(center);
              closest_vertex_p = vertex_p;
            }
          }

          if(std::abs(closest_vertex_p[1])<D && std::abs(closest_vertex_p[0])<D)
            criteria[ctr++] = 1;
          else
            criteria[ctr++] = 0;
        }
        GridRefinement::refine(coarse_grid, criteria, 0.5);
        coarse_grid.execute_coarsening_and_refinement();

        D /= 2.;
      }

      /*coarse_grid.refine_global(NUM_PRELIMINARY_GLOBAL_REF);

      for (unsigned int i = 0; i < NUM_PRELIMINARY_REF; ++i) {
        Vector<float> criteria(coarse_grid.n_active_cells());
        //cout  << "Active cells " << triangulation.n_active_cells() << endl;
        unsigned int ctr = 0;

        // Threshold
        const double max_thickness = 2. * l;
        const double min_thickness = 1.05 * l;
        const double D = min_thickness + (max_thickness-min_thickness)/(NUM_PRELIMINARY_REF-1)*(NUM_PRELIMINARY_REF-1-i);

        for (auto &cell : coarse_grid.active_cell_iterators()) {
          const Point<dim> c = cell->center();
          if(std::abs(c[1])<D && std::abs(c[0])<D)
            criteria[ctr++] = 1;
          else
            criteria[ctr++] = 0;
        }
        GridRefinement::refine(coarse_grid, criteria, 0.5);
        coarse_grid.execute_coarsening_and_refinement();
      }
      cout<<"Executed preliminary coarsening and refinement"<<endl;*/

    }

    // ------------------------------------------------------
    // Rectangle_1_99_manifold
    // ------------------------------------------------------

    class SquareManifold : public ChartManifold<2, 2>
    {
    public:
      SquareManifold(const Point<2> &center, const double half_length)
        : center(center), half_length(half_length) {}

      std::unique_ptr<Manifold<2>> clone() const override
      {
        return std::make_unique<SquareManifold>(center,half_length);
      }

      virtual Point<2> pull_back(const Point<2> &space_point) const override
      {
        // Map points to local coordinates within the square
        return Point<2>(space_point[0] - center[0], space_point[1] - center[1]);
      }

      virtual Point<2> push_forward(const Point<2> &chart_point) const override
      {
        // Ensure points are snapped to the square's boundary
        Point<2> result = center;

        result[0] += std::max(std::min(chart_point[0], half_length), -half_length);
        result[1] += std::max(std::min(chart_point[1], half_length), -half_length);
        return result;
      }

    private:
      const Point<2> center;
      const double half_length;
    };

    template <>
    void Rectangle_1_99_manifold<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      for (const auto &cell : coarse_grid.active_cell_iterators())
      {
        for (unsigned int face = 0; face < GeometryInfo<2>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary() && (cell->face(face)->boundary_id() == 1)) // Boundary ID 1 for the emitter, 9 for collector
          {
            cell->face(face)->set_manifold_id(1); // Assign manifold ID 1 for the emitter
          }
        }
      }

      const double l = 0.0001; // Square edge length
      const Point<2> center(0.0, 0.0); // Center of the square

      // Attach a square manifold to the emitter
      SquareManifold square_manifold(center, l );

      coarse_grid.set_manifold(1, square_manifold); // Set the manifold for the emitter
      coarse_grid.refine_global(3);

      /*const unsigned int NUM_PRELIMINARY_REF = 1;
      double D = 0.0033;
      for (unsigned int i = 0; i < NUM_PRELIMINARY_REF; ++i) {
        Vector<float> criteria(coarse_grid.n_active_cells());
        unsigned int ctr = 0;

        cout<<"D: "<<D<<endl;
        // Add some tolerance
        D *= 0.999;

        for (auto &cell : coarse_grid.active_cell_iterators()) {
          Point<2> closest_vertex_p = cell->center();
          double min_distance = closest_vertex_p.distance(center);
          for (const auto vertex : cell->vertex_indices()) {
            auto & vertex_p = cell->vertex(vertex);
            if (vertex_p.distance(center)<min_distance) {
              min_distance = vertex_p.distance(center);
              closest_vertex_p = vertex_p;
            }
          }

          if(std::abs(closest_vertex_p[1])<D && std::abs(closest_vertex_p[0])<D)
            criteria[ctr++] = 1;
          else
            criteria[ctr++] = 0;
        }
        GridRefinement::refine(coarse_grid, criteria, 0.5);
        coarse_grid.execute_coarsening_and_refinement();

        D /= 2.;
      }*/
    }


    // ------------------------------------------------------
    // angle_step14_forced
    // ------------------------------------------------------

    template <>
    void angle_step14_forced<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      coarse_grid.refine_global(1);
    }

    // ------------------------------------------------------
    // angle_Rectangle_1_100_forced
    // ------------------------------------------------------

    template <>
    void angle_Rectangle_1_100_forced<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);
    }

    // ------------------------------------------------------
    // Circular
    // ------------------------------------------------------

    template <>
    void Circular<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      const std::string path_to_mesh = "../mesh/cerchi_concentrici.msh";
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      /*double pi = 3.14159265358979323846;
      double Ve = 20000.;
      double l = 0.0004;
      double L = 0.004;
      cout<< "Exact flux: "<< 2 * pi * l * Ve / (L-l) <<std::endl;*/

    }

    // ------------------------------------------------------
    // LogCircular_1_2
    // ------------------------------------------------------

    template <>
    void LogCircular_1_2<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {

      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      const Point<2> center(0.0, 0.0);
      SphericalManifold<2> circular_manifold(center);


      coarse_grid.reset_all_manifolds();
      coarse_grid.set_all_manifold_ids(0);
      if (MANIFOLD_IS_APPLIED == 2 || MANIFOLD_IS_APPLIED == 3)
        coarse_grid.set_all_manifold_ids_on_boundary(1);
      else if (MANIFOLD_IS_APPLIED == 1)
        coarse_grid.set_all_manifold_ids(1);

      coarse_grid.set_manifold (1, circular_manifold);

      if (MANIFOLD_IS_APPLIED == 0 || MANIFOLD_IS_APPLIED == 2)
        coarse_grid.set_manifold (0, FlatManifold<2>());
      else if(MANIFOLD_IS_APPLIED == 3){
        TransfiniteInterpolationManifold<2> inner_manifold;
        inner_manifold.initialize(coarse_grid);
        coarse_grid.set_manifold (0, inner_manifold);
      }


      double l = 0.5;

      for (unsigned int i = 0; i < NUM_CONCENTRIC_REF; ++i) {
        Vector<float> criteria(coarse_grid.n_active_cells());
        //cout  << "Active cells " << triangulation.n_active_cells() << endl;
        unsigned int ctr = 0;

        // Threshold
        const double max_thickness = 1.5 * l;    // 2*
        const double min_thickness = 1.05 * l;    // 1.05
        double D = 0.;
        if (NUM_CONCENTRIC_REF==1)
          D = max_thickness;
        else
          D = min_thickness + (max_thickness-min_thickness)/(NUM_CONCENTRIC_REF-1)*(NUM_CONCENTRIC_REF-1-i);


        for (auto &cell : coarse_grid.active_cell_iterators()) {
          const Point<2> c = cell->center();
          if(center.distance(c) < D)
            criteria[ctr++] = 1;
          else
            criteria[ctr++] = 0;
        }
        GridRefinement::refine(coarse_grid, criteria, 0.5);
        coarse_grid.execute_coarsening_and_refinement();
        cout<<"Executed one concentric refinement"<<endl;
      }

      /*double pi = 3.14159265358979323846;
        double Ve = 20000.;
        double l = 0.5;
        double L = 1.0;
        cout<< std::scientific << std::setprecision(12)<< "Exact flux: "<< - (2 * pi * l) * (Ve / (log(l/L)*l)) <<std::endl;*/
    }

    // ------------------------------------------------------
    // LogCircular_1_10
    // ------------------------------------------------------

    template <>
    void LogCircular_1_10<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      //const std::string path_to_mesh = "../mesh/cerchi_concentrici_1_100.msh";
      //const std::string path_to_mesh = "../mesh/cerchi_concentrici.msh";
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      /*double pi = 3.14159265358979323846;
      double Ve = 20000.;
      double l = 0.0004;
      double L = 0.004;
      cout<< "Exact flux: "<< - (2 * pi * l) * (Ve / (log(l/L)*l)) <<std::endl;*/

      //ExactSolution exact_solution;
      //cout<<"ExactSolution at (0.0019, 0) = "<<exact_solution.value(Point<2>(0.0019, 0.),0)<< std::endl;

      const Point<2> center(0.0, 0.0);
      SphericalManifold<2> circular_manifold(center);
      coarse_grid.reset_all_manifolds();

      if (MANIFOLD_IS_APPLIED == 2)
        coarse_grid.set_all_manifold_ids_on_boundary(0);
      else if (MANIFOLD_IS_APPLIED == 1)
        coarse_grid.set_all_manifold_ids(0);

      coarse_grid.set_manifold (0, circular_manifold);


      double l = 0.0004;

      for (unsigned int i = 0; i < NUM_CONCENTRIC_REF; ++i) {

        // Threshold
        const double max_thickness = 1.5 * l;    // 2*
        const double min_thickness = 1.05 * l;    // 1.05
        double D = 0.;
        if (NUM_CONCENTRIC_REF==1)
          D = max_thickness;
        else
          D = min_thickness + (max_thickness-min_thickness)/(NUM_CONCENTRIC_REF-1)*(NUM_CONCENTRIC_REF-1-i);

        Vector<float> criteria(coarse_grid.n_active_cells());
        unsigned int ctr = 0;
        for (auto &cell : coarse_grid.active_cell_iterators()) {
          const Point<2> c = cell->center();
          if(center.distance(c) < D)
            criteria[ctr++] = 1;
          else
            criteria[ctr++] = 0;
        }
        GridRefinement::refine(coarse_grid, criteria, 0.5);
        coarse_grid.execute_coarsening_and_refinement();
        cout<<"Executed one concentric refinement"<<endl;
      }


    }

    // ------------------------------------------------------
    // LogCircular_1_100
    // ------------------------------------------------------

    template <>
    void LogCircular_1_100<2>::create_coarse_grid(Triangulation<2> &coarse_grid) {
      //const std::string path_to_mesh = "../mesh/cerchi_concentrici_1_100.msh";
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      /*double pi = 3.14159265358979323846;
      double Ve = 20000.;
      double l = 0.0004;
      double L = 0.04;
      cout << std::scientific << std::setprecision(12)<< "Exact flux: "<< - (2 * pi * l) * (Ve / (log(l/L)*l)) <<std::endl;*/

      //ExactSolution exact_solution;
      //cout<<"ExactSolution at (0.0019, 0) = "<<exact_solution.value(Point<2>(0.0019, 0.),0)<< std::endl;


      const Point<2> center(0.0, 0.0);
      SphericalManifold<2> circular_manifold(center);


      coarse_grid.reset_all_manifolds();
      coarse_grid.set_all_manifold_ids(0);
      if (MANIFOLD_IS_APPLIED == 2 || MANIFOLD_IS_APPLIED == 3)
        coarse_grid.set_all_manifold_ids_on_boundary(1);
      else if (MANIFOLD_IS_APPLIED == 1)
        coarse_grid.set_all_manifold_ids(1);

      coarse_grid.set_manifold (1, circular_manifold);

      if (MANIFOLD_IS_APPLIED == 2)
         coarse_grid.set_manifold (0, FlatManifold<2>());
      else if(MANIFOLD_IS_APPLIED == 3){
          TransfiniteInterpolationManifold<2> inner_manifold;
          inner_manifold.initialize(coarse_grid);
          coarse_grid.set_manifold (0, inner_manifold);
      }

    }




    // ------------------------------------------------------
    // CircularZeroDirichlet
    // ------------------------------------------------------

    template <>
    void CircularZeroDirichlet<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      //const std::string path_to_mesh = "../mesh/cerchi_concentrici_1_100.msh";
      const std::string path_to_mesh = "../mesh/cerchi_concentrici.msh";
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      /*double pi = 3.14159265358979323846;
      double l = 0.0004;
      double L = 0.004;
      std::cout << std::scientific << std::setprecision(12)
                << "Exact flux: " << - (2 * pi * l) * (pi / (L-l)) << std::endl;*/



      // Set up the circular manifold for the emitter (inner circle)
      const Point<2> center(0.0, 0.0); // Center of the circles

      for (const auto &cell : coarse_grid.active_cell_iterators())
      {
        for (unsigned int face = 0; face < GeometryInfo<2>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary() && (cell->face(face)->boundary_id() == 1 || cell->face(face)->boundary_id() == 9)) // Boundary ID 1 for the emitter, 9 for collector
          {
            cell->face(face)->set_manifold_id(1); // Assign manifold ID 1 for the emitter
          }
        }
      }

      // Attach a circular manifold to the emitter
      SphericalManifold<2> circular_manifold(center);
      coarse_grid.set_manifold(1, circular_manifold); // Set the manifold for the emitter

    }

    // ------------------------------------------------------
    // CircularStep14
    // ------------------------------------------------------

    template <>
    void CircularStep14<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      if (MANIFOLD_IS_APPLIED > 0){
        const Point<2> center(0.0, 0.0); // Center of the circles
        for (const auto &cell : coarse_grid.active_cell_iterators())
        {
          for (unsigned int face = 0; face < GeometryInfo<2>::faces_per_cell; ++face)
          {
            if (MANIFOLD_IS_APPLIED == 2) {
              if (cell->face(face)->at_boundary() && (cell->face(face)->boundary_id() == 1 || cell->face(face)->boundary_id() == 2 ||cell->face(face)->boundary_id() == 9)) // Boundary ID 1 for the emitter, 9 for collector
                cell->face(face)->set_manifold_id(1);
            }
            else if (MANIFOLD_IS_APPLIED == 1)
              cell->face(face)->set_manifold_id(1);
          }
        }

        // Attach a circular manifold to the emitter
        SphericalManifold<2> circular_manifold(center);
        coarse_grid.set_manifold(1, circular_manifold); // Set the manifold for the emitter
        }

    }

    template <int dim>
    void WireWire<dim>::create_coarse_grid(Triangulation<dim> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      GridOut grid_out;
      GridOutFlags::Msh msh_flags(true, true);
      grid_out.set_flags(msh_flags);
      grid_out.write_msh(coarse_grid, OUTPUT_PATH+"/initial_mesh.msh");

     const types::manifold_id emitterUp = 1;
      Point<2> center_emitterUp(0.0, 0.03);
      SphericalManifold<2> emitterUp_manifold(center_emitterUp);

      const types::manifold_id emitterDown = 2;
      Point<2> center_emitterDown(0.0, -0.03);
      SphericalManifold<2> emitterDown_manifold(center_emitterDown);

      const types::manifold_id collectorUp = 3;
      Point<2> center_collectorUp(0.053175, 0.03);
      SphericalManifold<2> collectorUp_manifold(center_collectorUp);

      const types::manifold_id collectorDown = 4;
      Point<2> center_collectorDown(0.053175, -0.03);
      SphericalManifold<2> collectorDown_manifold(center_collectorDown);

      TransfiniteInterpolationManifold<dim> inner_manifold;

      if (MANIFOLD_IS_APPLIED==1) {

        coarse_grid.reset_all_manifolds();
        coarse_grid.set_all_manifold_ids(0);

        const double R = 0.011;
        for (const auto &cell : coarse_grid.active_cell_iterators())
        {
          for (const auto &face : cell->face_iterators())
          {
            if (face->center().distance(center_emitterUp) < R)
              face->set_all_manifold_ids(emitterUp);
            else if (face->center().distance(center_emitterDown) < R)
              face->set_all_manifold_ids(emitterDown);
            else if (face->center().distance(center_collectorUp) < R)
              face->set_all_manifold_ids(collectorUp);
            else if (face->center().distance(center_collectorDown) < R)
              face->set_all_manifold_ids(collectorDown);

          }
        }
      } else if (MANIFOLD_IS_APPLIED==3) {
        coarse_grid.reset_all_manifolds();
        coarse_grid.set_all_manifold_ids(0);
        coarse_grid.set_all_manifold_ids_on_boundary(1, emitterUp);
        coarse_grid.set_all_manifold_ids_on_boundary(2, emitterDown);
        coarse_grid.set_all_manifold_ids_on_boundary(3, collectorUp);
        coarse_grid.set_all_manifold_ids_on_boundary(4, collectorDown);
      }

      coarse_grid.set_manifold(emitterUp, emitterUp_manifold);
      coarse_grid.set_manifold(emitterDown, emitterDown_manifold);
      coarse_grid.set_manifold(collectorUp, collectorUp_manifold);
      coarse_grid.set_manifold(collectorDown, collectorDown_manifold);

      inner_manifold.initialize(coarse_grid);
      coarse_grid.set_manifold (0, inner_manifold);
      //coarse_grid.set_manifold (0, FlatManifold<dim>());


      /*for (unsigned int i = 0; i < 2; ++i) {
        Vector<float> criteria(coarse_grid.n_active_cells());
        unsigned int ctr = 0;

        for (auto &cell : coarse_grid.active_cell_iterators()) {
          const Point<dim> c = cell->center();
          if (cell->center().distance(center_emitterUp) < R)
            criteria[ctr++] = 1;
          else if (cell->center().distance(center_emitterDown) < R)
            criteria[ctr++] = 1;
          else if (cell->center().distance(center_collectorUp) < R)
            criteria[ctr++] = 1;
          else if (cell->center().distance(center_collectorDown) < R)
            criteria[ctr++] = 1;
          else
            criteria[ctr++] = 0;
        }
        GridRefinement::refine(coarse_grid, criteria, 0.5);
        coarse_grid.execute_coarsening_and_refinement();
      }*/

      //coarse_grid.refine_global(1);
    }


    // Template instantiation
    template struct SetUpBase<2>;

    template struct CurvedRidges<2>;
    template struct SetUp<IonPropulsion::Data::CurvedRidges<2>, 2>;

    template struct SetupNone<2>;
    template struct SetUp<IonPropulsion::Data::SetupNone<2>, 2>;

    template struct Exercise_2_3<2>;
    template struct SetUp<IonPropulsion::Data::Exercise_2_3<2>, 2>;

    template struct Rectangle_1_99<2>;
    template struct SetUp<IonPropulsion::Data::Rectangle_1_99<2>, 2>;

    template struct Rectangle_1_99_manifold<2>;
    template struct SetUp<IonPropulsion::Data::Rectangle_1_99_manifold<2>, 2>;

    template struct angle_step14_forced<2>;
    template struct SetUp<IonPropulsion::Data::angle_step14_forced<2>, 2>;

    template struct angle_Rectangle_1_100_forced<2>;
    template struct SetUp<IonPropulsion::Data::angle_Rectangle_1_100_forced<2>, 2>;

    template struct Circular<2>;
    template struct SetUp<IonPropulsion::Data::Circular<2>, 2>;

    template struct LogCircular_1_2<2>;
    template struct SetUp<IonPropulsion::Data::LogCircular_1_2<2>, 2>;

    template struct LogCircular_1_10<2>;
    template struct SetUp<IonPropulsion::Data::LogCircular_1_10<2>, 2>;

    template struct LogCircular_1_100<2>;
    template struct SetUp<IonPropulsion::Data::LogCircular_1_100<2>, 2>;

    template struct CircularZeroDirichlet<2>;
    template struct SetUp<IonPropulsion::Data::CircularZeroDirichlet<2>, 2>;

    template struct CircularStep14<2>;
    template struct SetUp<IonPropulsion::Data::CircularStep14<2>, 2>;

    template struct WireWire<2>;
    template struct SetUp<IonPropulsion::Data::WireWire<2>, 2>;

  } // namespace Data
} // namespace IonPropulsion