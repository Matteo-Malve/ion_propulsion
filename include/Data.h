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

#ifndef DATA_H
#define DATA_H

#include "includes.h"
#include "CSVLogger.h"
#include "Constants.h"
#include <deal.II/grid/grid_in.h>


namespace IonPropulsion{
	using namespace dealii;

  namespace Data
  {
    constexpr double pi = 3.14159265358979323846;

    // ------------------------------------------------------
    // SetUpBase
    // ------------------------------------------------------
    template <int dim>
    struct SetUpBase : public Subscriptor
    {
      virtual const Function<dim> &get_boundary_values() const = 0;

      virtual const Function<dim> &get_right_hand_side() const = 0;

      virtual const Function<dim> &get_exact_solution() const = 0;

      virtual void
      create_coarse_grid( parallel::distributed::Triangulation<dim> &coarse_grid) const = 0;
    };

    // ------------------------------------------------------
    // SetUp
    // ------------------------------------------------------
    template <class Traits, int dim>
    struct SetUp : public SetUpBase<dim>
    {

      virtual const Function<dim> &get_boundary_values() const override;

      virtual const Function<dim> &get_right_hand_side() const override;

      virtual const Function<dim> &get_exact_solution() const override;

      virtual void
      create_coarse_grid( parallel::distributed::Triangulation<dim> &coarse_grid) const override;

    private:
      static const typename Traits::BoundaryValues boundary_values;
      static const typename Traits::RightHandSide  right_hand_side;
      static const typename Traits::ExactSolution  exact_solution;
    };

    // ------------------------------------------------------
    // Traits
    // ------------------------------------------------------

    template <class Traits, int dim>
    const typename Traits::BoundaryValues SetUp<Traits, dim>::boundary_values;
    template <class Traits, int dim>
    const typename Traits::RightHandSide SetUp<Traits, dim>::right_hand_side = typename Traits::RightHandSide();
    template <class Traits, int dim>
    const typename Traits::ExactSolution SetUp<Traits, dim>::exact_solution;

    // ------------------------------------------------------
    // SetupNone
    // ------------------------------------------------------

    template <int dim>
    struct SetupNone
    {
      class BoundaryValues : public Functions::ConstantFunction<dim>
      {
      public:
        BoundaryValues()
          : Functions::ConstantFunction<dim>(0.)
        {}

      };

      class RightHandSide : public Functions::ConstantFunction<dim>
      {
      public:
        RightHandSide()
          : Functions::ConstantFunction<dim>(1.)
        {}
      };

      class ExactSolution : public Functions::ConstantFunction<dim>
      {
      public:
        ExactSolution()
          : Functions::ConstantFunction<dim>(-1.e-19)  // Not available
        {}
      };

      static void create_coarse_grid( parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // Trait: DealiiStep14
    // ------------------------------------------------------

    template <int dim>
    struct DealiiStep14
    {
      // We need a class to denote the boundary values of the problem. In this
      // case, this is simple: it's the zero function, so don't even declare a
      // class, just an alias:
      //using BoundaryValues = Functions::ZeroFunction<dim>;
      class BoundaryValues : public Functions::ConstantFunction<dim>
      {
      public:
        BoundaryValues()
          : Functions::ConstantFunction<dim>(0.)
        {}

      };
      // Second, a class that denotes the right hand side. Since they are
      // constant, just subclass the corresponding class of the library and be
      // done:
      class RightHandSide : public Functions::ConstantFunction<dim>
      {
      public:
        RightHandSide()
          : Functions::ConstantFunction<dim>(1.)
        {}
      };

      class ExactSolution : public Functions::ConstantFunction<dim>
      {
      public:
        ExactSolution()
          : Functions::ConstantFunction<dim>(-1.e-19)  // Not available
        {}
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid( parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // LogCircular_1_2
    // ------------------------------------------------------

    template <int dim>
    struct LogCircular_1_2
    {

      class BoundaryValues : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;

          double Ve = 20000.;
          double l = 0.5;
          double L = 1.0;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return Ve / log(l/L) * log(r) - Ve * log(L) / log(l/L);
        }
      };

      class ExactSolution : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.5;
          double L = 1.0;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return Ve / log(l/L) * log(r) - Ve * log(L) / log(l/L);
        }
        virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.5;
          double L = 1.0;
          const auto x = p[0];
          const auto y = p[1];
          const double r2 = x*x + y*y;

          Tensor<1, dim> grad;
          grad[0] = Ve * x /(r2 * log(l/L));
          grad[1] = Ve * y /(r2 * log(l/L));
          return grad;
        };
      };

      using RightHandSide = Functions::ZeroFunction<dim>;

      static void create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid);
    };


    // ------------------------------------------------------
    // LogCircular_1_100
    // ------------------------------------------------------

    template <int dim>
    struct LogCircular_1_100
    {
      // We need a class to denote the boundary values of the problem. In this
      // case, this is simple: it's the zero function, so don't even declare a
      // class, just an alias:
      //using BoundaryValues = ExactSolution5b<dim>;
      class BoundaryValues : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.0004;
          // double L = 0.04;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return r < 1.1 * l ? Ve : 0.;
        }
      };
      class ExactSolution : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.0004;
          double L = 0.04;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return Ve / log(l/L) * log(r) - Ve * log(L) / log(l/L);
        }
        virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.0004;
          double L = 0.04;
          const auto x = p[0];
          const auto y = p[1];
          const double r2 = x*x + y*y;

          Tensor<1, dim> grad;
          grad[0] = Ve * x /(r2 * log(l/L));
          grad[1] = Ve * y /(r2 * log(l/L));
          return grad;
        };
      };

      using RightHandSide = Functions::ZeroFunction<dim>;

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // WireWire
    // ------------------------------------------------------

    template <int dim>
    struct WireWire
    {
      class BoundaryValues : public Functions::ConstantFunction<dim>
      {
      public:
        BoundaryValues()
          : Functions::ConstantFunction<dim>(0.)
        {}

      };

      class RightHandSide : public Functions::ConstantFunction<dim>
      {
      public:
        RightHandSide()
          : Functions::ConstantFunction<dim>(1.)
        {}
      };

      class ExactSolution : public Functions::ConstantFunction<dim>
      {
      public:
        ExactSolution()
          : Functions::ConstantFunction<dim>(-1.e-19)  // Not available
        {}
      };

      static void create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid);
    };


  } // namespace Data
}



#endif