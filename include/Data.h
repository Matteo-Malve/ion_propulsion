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
      create_coarse_grid(Triangulation<dim> &coarse_grid) const = 0;
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
      create_coarse_grid(Triangulation<dim> &coarse_grid) const override;

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
    const typename Traits::RightHandSide SetUp<Traits, dim>::right_hand_side;
    template <class Traits, int dim>
    const typename Traits::ExactSolution SetUp<Traits, dim>::exact_solution;


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

      static void create_coarse_grid(Triangulation<dim> &coarse_grid);
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
      static void create_coarse_grid(Triangulation<dim> &coarse_grid);
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
      static void create_coarse_grid(Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // FullTestSqruareComparison
    // ------------------------------------------------------
    constexpr double pi = 3.14159265358979323846;

    template <int dim>
    struct FullTestSqruareComparison
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

          double sigma2 = 0.0000005;
          double Ve = 20000.;
          double l = 0.0004;


          const auto x = p[0];
          const auto y = p[1];
          double r = sqrt(x*x + y*y);
          double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
          double freq = 0.5;
          double argX = (std::abs(x)-l) * pi / l * freq;
          double argY = (std::abs(y)-l) * pi / l * freq;

          return expTerm * Ve * (1 - sin(argX) * sin(argY));

        }
      };

      class ExactSolution : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;

          double sigma2 = 0.0000005;
          double Ve = 20000.;
          double l = 0.0004;


          const auto x = p[0];
          const auto y = p[1];
          double r = sqrt(x*x + y*y);
          double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
          double freq = 0.5;
          double argX = (std::abs(x)-l) * pi / l * freq;
          double argY = (std::abs(y)-l) * pi / l * freq;

          return expTerm * Ve * (1 - sin(argX) * sin(argY));

        }
      };

      class RightHandSide : public Function<dim>
      {
      public:
        virtual double value(const Point<dim> & p,
                             const unsigned int component) const override {
          (void)component;
          double sigma2 = 0.0000005;
          double freq = 0.5;
          double l = 0.0004;
          double Ve = 20000.;

          const auto x = p[0];
          const auto y = p[1];
          double r = sqrt(x*x + y*y);
          double expTerm = std::exp(- std::pow(r-l, 2) / (2 * sigma2));
          double argX = (std::abs(x)-l) * pi / l * freq;
          double argY = (std::abs(y)-l) * pi / l * freq;
          double signX = x>=0 ? +1. : -1.;
          double signY = y>=0 ? +1. : -1.;

          return - 1 * 1 *
                  (1 / ( l*l * sigma2*sigma2 * r)) * expTerm * Ve *
                  (   2 * freq * l * pi * sigma2 * x * (r-l) * cos(argX) * sin(argY) * signX +
                      freq*freq * pi*pi * sigma2*sigma2 * r * sin(argX) * sin(argY) +
                      2 * freq * l * pi * sigma2 * y * (r-l) * cos(argY) * sin(argX) * signY +
                      freq*freq * pi*pi * sigma2*sigma2 * r * sin(argX) * sin(argY) -
                      l * ( l * ( l*l * r + r * (r*r - 2*sigma2) +l * (sigma2 - 2*r*r) ) *  (- 1 + sin(argX) * sin(argY) ))
                  );
        }
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(Triangulation<dim> &coarse_grid);
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
      static void create_coarse_grid(Triangulation<dim> &coarse_grid);
    };

    // ------------------------------------------------------
    // LogCircular
    // ------------------------------------------------------

    template <int dim>
    struct LogCircular
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
          //double L = 0.04;

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
      };

      using RightHandSide = Functions::ZeroFunction<dim>;

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(Triangulation<dim> &coarse_grid);
    };


  } // namespace Data
}



#endif