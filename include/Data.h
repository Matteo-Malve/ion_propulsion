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
    // Trait: CurvedRidges
    // ------------------------------------------------------

    template <int dim>
    struct CurvedRidges
    {
      class BoundaryValues : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override;
      };


      class RightHandSide : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override;
      };

      class ExactSolution : public Functions::ConstantFunction<dim> {
      public:
        ExactSolution()
          : Functions::ConstantFunction<dim>(-1.e-19)  // Not available
        {}
      };

      static void create_coarse_grid( parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // Trait: Exercise_2_3
    // ------------------------------------------------------

    template <int dim>
    struct Exercise_2_3
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
    // Trait: Rectangle 1:99
    // ------------------------------------------------------

    template <int dim>
    struct Rectangle_1_99
    {
      // We need a class to denote the boundary values of the problem. In this
      // case, this is simple: it's the zero function, so don't even declare a
      // class, just an alias:
      using BoundaryValues = Functions::ZeroFunction<dim>;

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

      class ExactSolution : public Functions::ConstantFunction<dim> {
      public:
        ExactSolution()
          : Functions::ConstantFunction<dim>(-1.e-19)   // Not available
        {}
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid( parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // Trait: Rectangle_1_99_manifold
    // ------------------------------------------------------

    template <int dim>
    struct Rectangle_1_99_manifold
    {
      // We need a class to denote the boundary values of the problem. In this
      // case, this is simple: it's the zero function, so don't even declare a
      // class, just an alias:
      using BoundaryValues = Functions::ZeroFunction<dim>;

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

      class ExactSolution : public Functions::ConstantFunction<dim> {
      public:
        ExactSolution()
          : Functions::ConstantFunction<dim>(-1.e-19)   // Not available
        {}
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid( parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // angle_step14_forced
    // ------------------------------------------------------
    constexpr double pi = 3.14159265358979323846;

    template <int dim>
    struct angle_step14_forced
    {
      // We need a class to denote the boundary values of the problem. In this
      // case, this is simple: it's the zero function, so don't even declare a
      // class, just an alias:
      //using BoundaryValues = ExactSolution5b<dim>;
      class BoundaryValues : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double l = 0.5;
          double L = 1.0;

          const auto x = p[0];
          const auto y = p[1];

          return -(x-l)*(x-L)*(y-l)*(y-L)/(l*L)*x*y;
          //return -(x*y-l)*(x*y-L);
        }
      };

      class ExactSolution : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;

          double l = 0.5;
          double L = 1.0;

          const auto x = p[0];
          const auto y = p[1];

          return -(x-l)*(x-L)*(y-l)*(y-L)/(l*L)*x*y;
        }

        virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
          (void)component;

          double l = 0.5;
          double L = 1.0;
          const auto x = p[0];
          const auto y = p[1];

          Tensor<1, dim> grad;
          grad[0] = -1./(l*L) *
            ( l*(L-2.*x) + x*(3.*x-2.*L) )
            * (l-y) *(L-y) * y;
          grad[1] = 1./(l*L)*
            ( l*(L-2.*y) + y*(3.*y-2.*L) )
            * (l-x) *(x-L) * x;
          return grad;
        };
      };

      class RightHandSide : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double l = 0.5;
          double L = 1.0;
          const auto x = p[0];
          const auto y = p[1];
          double r2 = x*x+y*y;
          /*double signX = x>=0 ? +1. : -1.;
          double signY = y>=0 ? +1. : -1.;*/

          return - eps_0 * eps_r * 2./(l*L)*(
            L*(x+y)*(x+y)*(x+y) - L*L*r2 -3*x*y*r2
            + l*l * (L*(x+y)-r2)
            + l * ( L*L*(x+y) + (x+y)*(x+y)*(x+y) - 2*L*(x*x + 3*x*y + y*y ) )
          );


        }
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // angle_Rectangle_1_100_forced
    // ------------------------------------------------------
    template <int dim>
    struct angle_Rectangle_1_100_forced
    {

      class BoundaryValues : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double l = 0.0001;
          double L = 0.01;

          const auto x = p[0];
          const auto y = p[1];

          return -(x-l)*(x-L)*(y-l)*(y-L)/(l*L)*x*y;
          //return -(x*y-l)*(x*y-L);
        }
      };

      class ExactSolution : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;

          double l = 0.0001;
          double L = 0.01;

          const auto x = p[0];
          const auto y = p[1];

          return -(x-l)*(x-L)*(y-l)*(y-L)/(l*L)*x*y;
        }

        virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
          (void)component;

          double l = 0.0001;
          double L = 0.01;
          const auto x = p[0];
          const auto y = p[1];

          Tensor<1, dim> grad;
          grad[0] = -1./(l*L) *
            ( l*(L-2.*x) + x*(3.*x-2.*L) )
            * (l-y) *(L-y) * y;
          grad[1] = 1./(l*L)*
            ( l*(L-2.*y) + y*(3.*y-2.*L) )
            * (l-x) *(x-L) * x;
          return grad;
        };
      };

      class RightHandSide : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double l = 0.0001;
          double L = 0.01;
          const auto x = p[0];
          const auto y = p[1];
          double r2 = x*x+y*y;
          /*double signX = x>=0 ? +1. : -1.;
          double signY = y>=0 ? +1. : -1.;*/

          return - eps_0 * eps_r * 2./(l*L)*(
            L*(x+y)*(x+y)*(x+y) - L*L*r2 -3*x*y*r2
            + l*l * (L*(x+y)-r2)
            + l * ( L*L*(x+y) + (x+y)*(x+y)*(x+y) - 2*L*(x*x + 3*x*y + y*y ) )
          );


        }
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // Circular
    // ------------------------------------------------------

    template <int dim>
    struct Circular
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
          double L = 0.004;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return Ve * (r-L) / (l-L);
        }
      };
      class ExactSolution : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.0004;
          double L = 0.004;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return Ve * (r-L) / (l-L);

        }
      };

      class RightHandSide : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.0004;
          double L = 0.004;
          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          double laplacian = Ve / ((l-L)*r);
          return - eps_0 * eps_r * laplacian;
        }
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid);
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
          /*double Ve = 20000.;
          double l = 0.5;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return r < 1.1 * l ? Ve : 0.;*/

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
    // LogCircular_1_10
    // ------------------------------------------------------

    template <int dim>
    struct LogCircular_1_10
    {

      class BoundaryValues : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.0004;
          // double L = 0.004;

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
          double L = 0.004;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          return Ve / log(l/L) * log(r) - Ve * log(L) / log(l/L);
        }
        virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
          (void)component;
          double Ve = 20000.;
          double l = 0.0004;
          double L = 0.004;
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
    // CircularZeroDirichlet (1:10)
    // ------------------------------------------------------

    template <int dim>
    struct CircularZeroDirichlet
    {

      class BoundaryValues : public Functions::ConstantFunction<dim>
      {
      public:
        BoundaryValues()
          : Functions::ConstantFunction<dim>(0.)
        {}

      };

      class ExactSolution : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;

          double l = 0.0004;
          double L = 0.004;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);
          double arg = pi*(r-l)/(L-l);

          return sin(arg);
        }
        virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
          (void)component;

          double l = 0.0004;
          double L = 0.004;
          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);
          double arg = pi*(r-l)/(L-l);

          Tensor<1, dim> grad;
          grad[0] = (pi * x * cos(-arg) ) / ((L-l) * r);
          grad[1] = (pi * y * cos(-arg) ) / ((L-l) * r);
          return grad;
        };
      };

      class RightHandSide : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;

          double l = 0.0004;
          double L = 0.004;

          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);
          double arg = pi*(r-l)/(L-l);

          return eps_0 * eps_r *    //TODO: not sure of eps
            (pi * x*x * cos(arg) / ((L-l)* r*r*r))
          + (pi * y*y * cos(arg) / ((L-l)* r*r*r))
          - (2 * pi * cos(arg) / ((L-l)*r))
          + (pi*pi * x*x * sin(arg) / ((L-l)*(L-l)*r*r) )
          + (pi*pi * y*y * sin(arg) / ((L-l)*(L-l)*r*r));
        };
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(parallel::distributed::Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // CircularStep14
    // ------------------------------------------------------

    template <int dim>
    struct CircularStep14
    {

      class BoundaryValues : public Functions::ConstantFunction<dim>
      {
      public:
        BoundaryValues()
          : Functions::ConstantFunction<dim>(0.)
        {}

      };

      class ExactSolution : public Function<dim> {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double l = 0.5;
          double L = 1.0;

          const auto x = p[0];
          const auto y = p[1];
          const double r2 = (x*x + y*y);
          const double r = std::sqrt(r2);

          double C1 = - (L*L - l*l) / (4*log(l/L));
          double C2 = L*L/4 + log(L) * (L*L - l*l) / (4*log(l/L));

          return -r2/4. + C1 * log(r) + C2;
        }
        virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int component = 0) const override {
          (void)component;

          double l = 0.5;
          double L = 1.0;
          const auto x = p[0];
          const auto y = p[1];
          const double r = std::sqrt(x*x + y*y);

          double C1 = - (L*L - l*l) / (4*log(l/L));

          double u_r = - r / 2. + C1 / r;

          Tensor<1, dim> grad;
          grad[0] = u_r * x / r;
          grad[1] = u_r * y / r;
          return grad;
        };
      };

      class RightHandSide : public Functions::ConstantFunction<dim>
      {
      public:
        RightHandSide()
          : Functions::ConstantFunction<dim>(1.)
        {}
      };

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