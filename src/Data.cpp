/* ------------------------------------------------------------------------
*
 * SPDX-License-Identifier: GPL-3-or-later
 * Copyright (C) 2023 - 2024 by Matteo Malvestiti
 *
 * This file is part of ion_propulsion.
 *
 * ion_propulsion is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ion_propulsion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ion_propulsion; see the file COPYING.  If not, see
 * <https://www.gnu.org/licenses/>.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Matteo Malvestiti, Politecnico di Milano, 2024
 *
 */

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
      parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      Traits::create_coarse_grid(coarse_grid);
    }

    // ------------------------------------------------------
    // SetupNone
    // ------------------------------------------------------

    template <int dim>
    void SetupNone<dim>::create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      // If needed, define your manifold description HERE

    }


    // ------------------------------------------------------
    // DealiiStep14
    // ------------------------------------------------------

    template <>
    void DealiiStep14<2>::create_coarse_grid(parallel::distributed::Triangulation<2> &coarse_grid)
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


      std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
      for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      }

      coarse_grid.create_triangulation(vertices, cells, SubCellData());

      for(auto cell : coarse_grid.active_cell_iterators())
        for (unsigned int face = 0; face < GeometryInfo<2>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(9);

      coarse_grid.refine_global(1);
    }

    // ------------------------------------------------------
    // LogCircular_1_2
    // ------------------------------------------------------

    template <>
    void LogCircular_1_2<2>::create_coarse_grid(parallel::distributed::Triangulation<2> &coarse_grid)
    {

      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

      const Point<2> center(0.0, 0.0); // Center of the circles

      if (MANIFOLD_IS_APPLIED>0){
        for (const auto &cell : coarse_grid.active_cell_iterators())
        {
          for (unsigned int face = 0; face < GeometryInfo<2>::faces_per_cell; ++face)
          {
            if (MANIFOLD_IS_APPLIED==2) {
              if (cell->face(face)->at_boundary() && (cell->face(face)->boundary_id() == 1 || cell->face(face)->boundary_id() == 3 ||cell->face(face)->boundary_id() == 9)) // Boundary ID 1 for the emitter, 9 for collector
              {
                cell->face(face)->set_manifold_id(1); // Assign manifold ID 1 for the emitter
              }
            } else if (MANIFOLD_IS_APPLIED == 1) {
              cell->face(face)->set_manifold_id(1);
            } else {
              DEAL_II_NOT_IMPLEMENTED();
            }
          }
        }

        // Attach a circular manifold to the emitter
        SphericalManifold<2> circular_manifold(center);
        coarse_grid.set_manifold(1, circular_manifold); // Set the manifold for the emitter
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
    // LogCircular_1_100
    // ------------------------------------------------------

    template <>
    void LogCircular_1_100<2>::create_coarse_grid(parallel::distributed::Triangulation<2> &coarse_grid) {
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
    // WireWire
    // ------------------------------------------------------

    template <int dim>
    void WireWire<dim>::create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid)
    {
      const std::string path_to_mesh = PATH_TO_MESH;
      std::ifstream input_file(path_to_mesh);
      GridIn<2>       grid_in;
      grid_in.attach_triangulation(coarse_grid);
      grid_in.read_msh(input_file);

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
      } else if (MANIFOLD_IS_APPLIED==2 || MANIFOLD_IS_APPLIED==3) {
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

      if (MANIFOLD_IS_APPLIED==3) {
        inner_manifold.initialize(coarse_grid);
        coarse_grid.set_manifold (0, inner_manifold);
      }
      if (MANIFOLD_IS_APPLIED==2)
        coarse_grid.set_manifold (0, FlatManifold<dim>());


    }


    // Template instantiation
    template struct SetUpBase<2>;

    template struct SetupNone<2>;
    template struct SetUp<IonPropulsion::Data::SetupNone<2>, 2>;

    template struct DealiiStep14<2>;
    template struct SetUp<IonPropulsion::Data::DealiiStep14<2>, 2>;

    template struct LogCircular_1_2<2>;
    template struct SetUp<IonPropulsion::Data::LogCircular_1_2<2>, 2>;

    template struct LogCircular_1_100<2>;
    template struct SetUp<IonPropulsion::Data::LogCircular_1_100<2>, 2>;

    template struct WireWire<2>;
    template struct SetUp<IonPropulsion::Data::WireWire<2>, 2>;

  } // namespace Data
} // namespace IonPropulsion