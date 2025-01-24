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
      virtual void update_convergence_table() = 0;
      virtual void print_convergence_table() const {};

      template <class Archive>
      void serialize(Archive &ar, const unsigned int version) {
        (void)ar;
        (void)version;
      };

      void checkpoint();
      void restart();

    protected:
      const SmartPointer<Triangulation<dim>> triangulation;
      std::shared_ptr<ConvergenceTable> convergence_table;

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
             const Function<dim> &      boundary_values,
             const unsigned degree
             );
      virtual ~Solver() override;

      virtual void solve_problem() override;

      virtual void postprocess(
        const Evaluation::EvaluationBase<dim> &postprocessor) const override;

      virtual unsigned int n_dofs() const override;

      void update_convergence_table() override;

    protected:
      const SmartPointer<const FiniteElement<dim>>  fe;
      const SmartPointer<const Quadrature<dim>>     quadrature;
      const SmartPointer<const Quadrature<dim - 1>> face_quadrature;
      DoFHandler<dim>                               dof_handler;
      Vector<double>                                solution;
      Vector<double>                                homogeneous_solution;
      Vector<double>                                Rg_vector;
      const SmartPointer<const Function<dim>>       boundary_values;
      double                                        conservative_flux;
      MappingQ<dim>      mapping;

      virtual void assemble_rhs(Vector<double> &rhs) const = 0;

      virtual void construct_Rg_vector() = 0;

      virtual void retrieve_Rg() = 0;


      struct LinearSystem
      {
        LinearSystem(const DoFHandler<dim> &dof_handler);

        void solve(Vector<double> &solution) const;

        AffineConstraints<double> hanging_node_constraints;
        SparsityPattern           sparsity_pattern;
        SparseMatrix<double>      matrix;
        Vector<double>            rhs;
        SparseMatrix<double>      Umatrix;
      };

      std::unique_ptr<LinearSystem>                 linear_system_ptr;

      virtual void compute_second_order_flux(LinearSystem &linear_system) = 0;

    private:

      struct AssemblyScratchData
      {
        AssemblyScratchData(const FiniteElement<dim> &fe,
                            const Quadrature<dim> &   quadrature,
                            MappingQ<dim> & mapping);
        AssemblyScratchData(const AssemblyScratchData &scratch_data);

        MappingQ<dim>  mapping;
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
                   const Function<dim> &      boundary_values,
                   const unsigned degree);

      virtual void output_solution() const override;

    protected:
      const SmartPointer<const Function<dim>>       rhs_function;
      double                                        conservative_flux;
      Vector<double>                                Au;

      virtual void assemble_rhs(Vector<double> &rhs) const override;

      virtual void construct_Rg_vector() override;

      void compute_second_order_flux(typename Solver<dim>::LinearSystem &linear_system) override;

    private:
      void retrieve_Rg() override {
        this->solution += this->Rg_vector;
      }
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
        const DualFunctional::DualFunctionalBase<dim> &dual_functional,
        const Function<dim> &                          special_rhs_function,
        const Function<dim> &                          special_boundary_values,
        const unsigned degree);

    protected:
      const SmartPointer<const Function<dim>>       special_rhs_function;
      const SmartPointer<const Function<dim>>       special_boundary_values;

      const SmartPointer<const DualFunctional::DualFunctionalBase<dim>>
                   dual_functional;
      virtual void assemble_rhs(Vector<double> &rhs) const override;  //TODO: do it in WeightedResidual
      virtual void conservative_flux_rhs(Vector<double> & rhs) const = 0;   //TODO: New
      static const Functions::ZeroFunction<dim> boundary_values;

      /*void assemble_conservative_flux_rhs(
        Vector<double> &rhs
        ,SparseMatrix<double> & Umatrix
        );*/

    private:
      virtual void construct_Rg_vector() override {};
      void retrieve_Rg() override {};

      void compute_second_order_flux(typename Solver<dim>::LinearSystem & ) override {};

    };

    template <int dim>
    const Functions::ZeroFunction<dim> DualSolver<dim>::boundary_values;



  } // namespace LaplaceSolver
}

#endif