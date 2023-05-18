#ifndef GETPOT_GOE_DATA_H
#define GETPOT_GOE_DATA_H

#include "GOE_LaplaceSolver.h"

namespace GOE{
    using namespace dealii;
    namespace Data
    {

        template <int dim>
        struct SetUpBase : public Subscriptor
        {
            virtual const Function<dim> &get_boundary_values() const = 0;

            virtual const Function<dim> &get_right_hand_side() const = 0;

            virtual void
            create_coarse_grid(Triangulation<dim> &coarse_grid) const = 0;
        };


        template <class Traits, int dim>
        struct SetUp : public SetUpBase<dim>
        {
            virtual const Function<dim> &get_boundary_values() const override;

            virtual const Function<dim> &get_right_hand_side() const override;


            virtual void
            create_coarse_grid(Triangulation<dim> &coarse_grid) const override;

        private:
            static const typename Traits::BoundaryValues boundary_values;
            static const typename Traits::RightHandSide  right_hand_side;
        };

        template <class Traits, int dim>
        const typename Traits::BoundaryValues SetUp<Traits, dim>::boundary_values;
        template <class Traits, int dim>
        const typename Traits::RightHandSide SetUp<Traits, dim>::right_hand_side;

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



        template <int dim>
        struct CurvedRidges
        {
            class BoundaryValues : public Function<dim>
            {
            public:
                virtual double value(const Point<dim> & p,
                                     const unsigned int component) const;
            };


            class RightHandSide : public Function<dim>
            {
            public:
                virtual double value(const Point<dim> & p,
                                     const unsigned int component) const;
            };

            static void create_coarse_grid(Triangulation<dim> &coarse_grid);
        };


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



        template <int dim>
        struct Exercise_2_3
        {
            using BoundaryValues = Functions::ZeroFunction<dim>;

            class RightHandSide : public Functions::ConstantFunction<dim>
            {
            public:
                RightHandSide()
                        : Functions::ConstantFunction<dim>(1.)
                {}
            };

            static void create_coarse_grid(Triangulation<dim> &coarse_grid);
        };


        template <>
        void Exercise_2_3<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
        {
            const unsigned int dim = 2;

            const std::vector<Point<2>> vertices = {
                    {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0},
                    {-1.0, -0.5}, {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5},
                    {-1.0, +0.0}, {-0.5, +0.0}, {+0.5, +0.0}, {+1.0, +0.0},
                    {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5}, {+1.0, +0.5},
                    {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}};

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

            coarse_grid.refine_global(1);
        }
    } // namespace Data
}

#endif //GETPOT_GOE_DATA_H
