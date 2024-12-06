#ifndef LAPLACE_SOLVER_H
#define LAPLACE_SOLVER_H

#include "includes.h"
#include "Evaluation.h"
#include "DualFunctional.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{

    // ------------------------------------------------------
    // Base
    // ------------------------------------------------------
    template <int dim>
    class Base
    {
    public:
      Base(Triangulation<dim> &coarse_grid);
      virtual ~Base() = default;

      virtual void solve_problem() = 0;
      virtual void postprocess(
        const Evaluation::EvaluationBase<dim> &postprocessor) const = 0;
      virtual void         refine_grid()                            = 0;
      virtual unsigned int n_dofs() const                           = 0;

      virtual void set_refinement_cycle(const unsigned int cycle);

      virtual void output_solution() const = 0;

    protected:
      const SmartPointer<Triangulation<dim>> triangulation;

      unsigned int refinement_cycle;
    };

    // ------------------------------------------------------
    // Solver
    // ------------------------------------------------------
    template <int dim>
    class Solver : public virtual Base<dim>
    {
    public:
      Solver(Triangulation<dim> &       triangulation,
             const FiniteElement<dim> & fe,
             const Quadrature<dim> &    quadrature,
             const Quadrature<dim - 1> &face_quadrature,
             const Function<dim> &      boundary_values);
      virtual ~Solver() override;

      virtual void solve_problem() override;

      virtual void postprocess(
        const Evaluation::EvaluationBase<dim> &postprocessor) const override;

      virtual unsigned int n_dofs() const override;

    protected:
      const SmartPointer<const FiniteElement<dim>>  fe;
      const SmartPointer<const Quadrature<dim>>     quadrature;
      const SmartPointer<const Quadrature<dim - 1>> face_quadrature;
      DoFHandler<dim>                               dof_handler;
      Vector<double>                                solution;
      const SmartPointer<const Function<dim>>       boundary_values;

      virtual void assemble_rhs(Vector<double> &rhs) const = 0;

    private:
      struct LinearSystem
      {
        LinearSystem(const DoFHandler<dim> &dof_handler);

        void solve(Vector<double> &solution) const;

        AffineConstraints<double> hanging_node_constraints;
        SparsityPattern           sparsity_pattern;
        SparseMatrix<double>      matrix;
        Vector<double>            rhs;
      };


      struct AssemblyScratchData
      {
        AssemblyScratchData(const FiniteElement<dim> &fe,
                            const Quadrature<dim> &   quadrature);
        AssemblyScratchData(const AssemblyScratchData &scratch_data);

        FEValues<dim> fe_values;
      };

      struct AssemblyCopyData
      {
        FullMatrix<double>                   cell_matrix;
        std::vector<types::global_dof_index> local_dof_indices;
      };


      void assemble_linear_system(LinearSystem &linear_system);

      void local_assemble_matrix(
        const typename DoFHandler<dim>::active_cell_iterator &cell,
        AssemblyScratchData &                                 scratch_data,
        AssemblyCopyData &                                    copy_data) const;


      void copy_local_to_global(const AssemblyCopyData &copy_data,
                                LinearSystem &          linear_system) const;
    };


    // ------------------------------------------------------
    // PrimalSolver
    // ------------------------------------------------------
    template <int dim>
    class PrimalSolver : public Solver<dim>
    {
    public:
      PrimalSolver(Triangulation<dim> &       triangulation,
                   const FiniteElement<dim> & fe,
                   const Quadrature<dim> &    quadrature,
                   const Quadrature<dim - 1> &face_quadrature,
                   const Function<dim> &      rhs_function,
                   const Function<dim> &      boundary_values);

      virtual void output_solution() const override;

    protected:
      const SmartPointer<const Function<dim>> rhs_function;
      virtual void assemble_rhs(Vector<double> &rhs) const override;
    };

    // ------------------------------------------------------
    // DualSolver
    // ------------------------------------------------------
    template <int dim>
    class DualSolver : public Solver<dim>
    {
    public:
      DualSolver(
        Triangulation<dim> &                           triangulation,
        const FiniteElement<dim> &                     fe,
        const Quadrature<dim> &                        quadrature,
        const Quadrature<dim - 1> &                    face_quadrature,
        const DualFunctional::DualFunctionalBase<dim> &dual_functional);

    protected:
      const SmartPointer<const DualFunctional::DualFunctionalBase<dim>>
                   dual_functional;
      virtual void assemble_rhs(Vector<double> &rhs) const override;

      static const Functions::ZeroFunction<dim> boundary_values;
    };

    template <int dim>
    const Functions::ZeroFunction<dim> DualSolver<dim>::boundary_values;



  } // namespace LaplaceSolver
}

#endif