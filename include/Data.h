#ifndef DATA_H
#define DATA_H

#include "includes.h"

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


      virtual void
      create_coarse_grid(Triangulation<dim> &coarse_grid) const override;

    private:
      static const typename Traits::BoundaryValues boundary_values;
      static const typename Traits::RightHandSide  right_hand_side;
    };

    // ------------------------------------------------------
    // Traits
    // ------------------------------------------------------

    template <class Traits, int dim>
    const typename Traits::BoundaryValues SetUp<Traits, dim>::boundary_values;
    template <class Traits, int dim>
    const typename Traits::RightHandSide SetUp<Traits, dim>::right_hand_side;


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
          : Functions::ConstantFunction<dim>(2.)
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

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static void create_coarse_grid(Triangulation<dim> &coarse_grid);
    };
  } // namespace Data
}



#endif