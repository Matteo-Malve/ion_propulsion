#ifndef LAPLACE_SOLVER_H
#define LAPLACE_SOLVER_H

#include "includes.h"

namespace IonPropulsion{
  using namespace dealii;
  namespace LaplaceSolver{
  // @sect4{The Laplace solver base class}

  // This class is almost unchanged, with the exception that it declares two
  // more functions: <code>output_solution</code> will be used to generate
  // output files from the actual solutions computed by derived classes, and
  // the <code>set_refinement_cycle</code> function by which the testing
  // framework sets the number of the refinement cycle to a local variable
  // in this class; this number is later used to generate filenames for the
  // solution output.
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


  template <int dim>
  Base<dim>::Base(Triangulation<dim> &coarse_grid)
    : triangulation(&coarse_grid)
    , refinement_cycle(numbers::invalid_unsigned_int)
  {}



  template <int dim>
  void Base<dim>::set_refinement_cycle(const unsigned int cycle)
  {
    refinement_cycle = cycle;
  }


  // @sect4{The Laplace Solver class}

  // Likewise, the <code>Solver</code> class is entirely unchanged and will
  // thus not be discussed.
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


    // The remainder of the class is essentially a copy of step-13
    // as well, including the data structures and functions
    // necessary to compute the linear system in parallel using the
    // WorkStream framework:
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



  template <int dim>
  Solver<dim>::Solver(Triangulation<dim> &       triangulation,
                      const FiniteElement<dim> & fe,
                      const Quadrature<dim> &    quadrature,
                      const Quadrature<dim - 1> &face_quadrature,
                      const Function<dim> &      boundary_values)
    : Base<dim>(triangulation)
    , fe(&fe)
    , quadrature(&quadrature)
    , face_quadrature(&face_quadrature)
    , dof_handler(triangulation)
    , boundary_values(&boundary_values)
  {}


  template <int dim>
  Solver<dim>::~Solver()
  {
    dof_handler.clear();
  }


  template <int dim>
  void Solver<dim>::solve_problem()
  {
    dof_handler.distribute_dofs(*fe);
    solution.reinit(dof_handler.n_dofs());

    LinearSystem linear_system(dof_handler);
    assemble_linear_system(linear_system);
    linear_system.solve(solution);
  }


  template <int dim>
  void Solver<dim>::postprocess(
    const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    postprocessor(dof_handler, solution);
  }


  template <int dim>
  unsigned int Solver<dim>::n_dofs() const
  {
    return dof_handler.n_dofs();
  }


  // The following few functions and constructors are verbatim
  // copies taken from step-13:
  template <int dim>
  void Solver<dim>::assemble_linear_system(LinearSystem &linear_system)
  {
    Threads::Task<void> rhs_task =
      Threads::new_task(&Solver<dim>::assemble_rhs, *this, linear_system.rhs);

    auto worker =
      [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
             AssemblyScratchData &scratch_data,
             AssemblyCopyData &   copy_data) {
        this->local_assemble_matrix(cell, scratch_data, copy_data);
      };

    auto copier = [this, &linear_system](const AssemblyCopyData &copy_data) {
      this->copy_local_to_global(copy_data, linear_system);
    };

    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    worker,
                    copier,
                    AssemblyScratchData(*fe, *quadrature),
                    AssemblyCopyData());
    linear_system.hanging_node_constraints.condense(linear_system.matrix);

    std::map<types::global_dof_index, double> boundary_value_map;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             *boundary_values,
                                             boundary_value_map);

    rhs_task.join();
    linear_system.hanging_node_constraints.condense(linear_system.rhs);

    MatrixTools::apply_boundary_values(boundary_value_map,
                                       linear_system.matrix,
                                       solution,
                                       linear_system.rhs);
  }


  template <int dim>
  Solver<dim>::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<dim> &fe,
    const Quadrature<dim> &   quadrature)
    : fe_values(fe, quadrature, update_gradients | update_JxW_values)
  {}


  template <int dim>
  Solver<dim>::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : fe_values(scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_gradients | update_JxW_values)
  {}


  template <int dim>
  void Solver<dim>::local_assemble_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    AssemblyScratchData &                                 scratch_data,
    AssemblyCopyData &                                    copy_data) const
  {
    const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature->size();

    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);

    copy_data.local_dof_indices.resize(dofs_per_cell);

    scratch_data.fe_values.reinit(cell);

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          copy_data.cell_matrix(i, j) +=
            (scratch_data.fe_values.shape_grad(i, q_point) *
             scratch_data.fe_values.shape_grad(j, q_point) *
             scratch_data.fe_values.JxW(q_point));

    cell->get_dof_indices(copy_data.local_dof_indices);
  }



  template <int dim>
  void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data,
                                         LinearSystem &linear_system) const
  {
    for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i)
      for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j)
        linear_system.matrix.add(copy_data.local_dof_indices[i],
                                 copy_data.local_dof_indices[j],
                                 copy_data.cell_matrix(i, j));
  }


  // Now for the functions that implement actions in the linear
  // system class. First, the constructor initializes all data
  // elements to their correct sizes, and sets up a number of
  // additional data structures, such as constraints due to hanging
  // nodes. Since setting up the hanging nodes and finding out about
  // the nonzero elements of the matrix is independent, we do that
  // in parallel (if the library was configured to use concurrency,
  // at least; otherwise, the actions are performed
  // sequentially). Note that we start only one thread, and do the
  // second action in the main thread. Since only one thread is
  // generated, we don't use the <code>Threads::TaskGroup</code>
  // class here, but rather use the one created task object
  // directly to wait for this particular task's exit. The
  // approach is generally the same as the one we have used in
  // <code>Solver::assemble_linear_system()</code> above.
  //
  // Note that taking the address of the
  // <code>DoFTools::make_hanging_node_constraints</code> function
  // is a little tricky, since there are actually three functions of
  // this name, one for each supported space dimension. Taking
  // addresses of overloaded functions is somewhat complicated in
  // C++, since the address-of operator <code>&</code> in that case
  // returns a set of values (the addresses of all
  // functions with that name), and selecting the right one is then
  // the next step. If the context dictates which one to take (for
  // example by assigning to a function pointer of known type), then
  // the compiler can do that by itself, but if this set of pointers
  // shall be given as the argument to a function that takes a
  // template, the compiler could choose all without having a
  // preference for one. We therefore have to make it clear to the
  // compiler which one we would like to have; for this, we could
  // use a cast, but for more clarity, we assign it to a temporary
  // <code>mhnc_p</code> (short for <code>pointer to
  // make_hanging_node_constraints</code>) with the right type, and
  // using this pointer instead.
  template <int dim>
  Solver<dim>::LinearSystem::LinearSystem(const DoFHandler<dim> &dof_handler)
  {
    hanging_node_constraints.clear();

    void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) =
      &DoFTools::make_hanging_node_constraints;

    // Start a side task then continue on the main thread
    Threads::Task<void> side_task =
      Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);



    // Wait for the side task to be done before going further
    side_task.join();

    hanging_node_constraints.close();
    hanging_node_constraints.condense(dsp);
    sparsity_pattern.copy_from(dsp);

    matrix.reinit(sparsity_pattern);
    rhs.reinit(dof_handler.n_dofs());
  }



  template <int dim>
  void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const
  {
    SolverControl            solver_control(5000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(matrix, 1.2);

    cg.solve(matrix, solution, rhs, preconditioner);

    hanging_node_constraints.distribute(solution);
  }



  // @sect4{The PrimalSolver class}

  // The <code>PrimalSolver</code> class is also mostly unchanged except for
  // implementing the <code>output_solution</code> function. We keep the
  // <code>GlobalRefinement</code> and <code>RefinementKelly</code> classes
  // in this program, and they can then rely on the default implementation
  // of this function which simply outputs the primal solution. The class
  // implementing dual weighted error estimators will overload this function
  // itself, to also output the dual solution.
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


  template <int dim>
  PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &       triangulation,
                                  const FiniteElement<dim> & fe,
                                  const Quadrature<dim> &    quadrature,
                                  const Quadrature<dim - 1> &face_quadrature,
                                  const Function<dim> &      rhs_function,
                                  const Function<dim> &      boundary_values)
    : Base<dim>(triangulation)
    , Solver<dim>(triangulation,
                  fe,
                  quadrature,
                  face_quadrature,
                  boundary_values)
    , rhs_function(&rhs_function)
  {}



  template <int dim>
  void PrimalSolver<dim>::output_solution() const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "solution");
    data_out.build_patches();

    std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                      ".vtu");
    data_out.write(out, DataOutBase::vtu);
  }



  template <int dim>
  void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const
  {
    FEValues<dim> fe_values(*this->fe,
                            *this->quadrature,
                            update_values | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
    const unsigned int n_q_points    = this->quadrature->size();

    Vector<double>                       cell_rhs(dofs_per_cell);
    std::vector<double>                  rhs_values(n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : this->dof_handler.active_cell_iterators())
      {
        cell_rhs = 0;

        fe_values.reinit(cell);

        rhs_function->value_list(fe_values.get_quadrature_points(),
                                 rhs_values);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                            rhs_values[q_point] *               // f((x_q)
                            fe_values.JxW(q_point));            // dx

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          rhs(local_dof_indices[i]) += cell_rhs(i);
      }
  }


  // @sect4{The RefinementGlobal and RefinementKelly classes}

  // For the following two classes, the same applies as for most of the
  // above: the class is taken from the previous example as-is:
  template <int dim>
  class RefinementGlobal : public PrimalSolver<dim>
  {
  public:
    RefinementGlobal(Triangulation<dim> &       coarse_grid,
                     const FiniteElement<dim> & fe,
                     const Quadrature<dim> &    quadrature,
                     const Quadrature<dim - 1> &face_quadrature,
                     const Function<dim> &      rhs_function,
                     const Function<dim> &      boundary_values);

    virtual void refine_grid() override;
  };



  template <int dim>
  RefinementGlobal<dim>::RefinementGlobal(
    Triangulation<dim> &       coarse_grid,
    const FiniteElement<dim> & fe,
    const Quadrature<dim> &    quadrature,
    const Quadrature<dim - 1> &face_quadrature,
    const Function<dim> &      rhs_function,
    const Function<dim> &      boundary_values)
    : Base<dim>(coarse_grid)
    , PrimalSolver<dim>(coarse_grid,
                        fe,
                        quadrature,
                        face_quadrature,
                        rhs_function,
                        boundary_values)
  {}



  template <int dim>
  void RefinementGlobal<dim>::refine_grid()
  {
    this->triangulation->refine_global(1);
  }



  template <int dim>
  class RefinementKelly : public PrimalSolver<dim>
  {
  public:
    RefinementKelly(Triangulation<dim> &       coarse_grid,
                    const FiniteElement<dim> & fe,
                    const Quadrature<dim> &    quadrature,
                    const Quadrature<dim - 1> &face_quadrature,
                    const Function<dim> &      rhs_function,
                    const Function<dim> &      boundary_values);

    virtual void refine_grid() override;
  };



  template <int dim>
  RefinementKelly<dim>::RefinementKelly(
    Triangulation<dim> &       coarse_grid,
    const FiniteElement<dim> & fe,
    const Quadrature<dim> &    quadrature,
    const Quadrature<dim - 1> &face_quadrature,
    const Function<dim> &      rhs_function,
    const Function<dim> &      boundary_values)
    : Base<dim>(coarse_grid)
    , PrimalSolver<dim>(coarse_grid,
                        fe,
                        quadrature,
                        face_quadrature,
                        rhs_function,
                        boundary_values)
  {}



  template <int dim>
  void RefinementKelly<dim>::refine_grid()
  {
    Vector<float> estimated_error_per_cell(
      this->triangulation->n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      this->dof_handler,
      QGauss<dim - 1>(this->fe->degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      this->solution,
      estimated_error_per_cell);
    GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);
    this->triangulation->execute_coarsening_and_refinement();
  }



  // @sect4{The RefinementWeightedKelly class}

  // This class is a variant of the previous one, in that it allows to
  // weight the refinement indicators we get from the library's Kelly
  // indicator by some function. We include this class since the goal of
  // this example program is to demonstrate automatic refinement criteria
  // even for complex output quantities such as point values or stresses. If
  // we did not solve a dual problem and compute the weights thereof, we
  // would probably be tempted to give a hand-crafted weighting to the
  // indicators to account for the fact that we are going to evaluate these
  // quantities. This class accepts such a weighting function as argument to
  // its constructor:
  template <int dim>
  class RefinementWeightedKelly : public PrimalSolver<dim>
  {
  public:
    RefinementWeightedKelly(Triangulation<dim> &       coarse_grid,
                            const FiniteElement<dim> & fe,
                            const Quadrature<dim> &    quadrature,
                            const Quadrature<dim - 1> &face_quadrature,
                            const Function<dim> &      rhs_function,
                            const Function<dim> &      boundary_values,
                            const Function<dim> &      weighting_function);

    virtual void refine_grid() override;

  private:
    const SmartPointer<const Function<dim>> weighting_function;
  };



  template <int dim>
  RefinementWeightedKelly<dim>::RefinementWeightedKelly(
    Triangulation<dim> &       coarse_grid,
    const FiniteElement<dim> & fe,
    const Quadrature<dim> &    quadrature,
    const Quadrature<dim - 1> &face_quadrature,
    const Function<dim> &      rhs_function,
    const Function<dim> &      boundary_values,
    const Function<dim> &      weighting_function)
    : Base<dim>(coarse_grid)
    , PrimalSolver<dim>(coarse_grid,
                        fe,
                        quadrature,
                        face_quadrature,
                        rhs_function,
                        boundary_values)
    , weighting_function(&weighting_function)
  {}



  // Now, here comes the main function, including the weighting:
  template <int dim>
  void RefinementWeightedKelly<dim>::refine_grid()
  {
    // First compute some residual based error indicators for all cells by a
    // method already implemented in the library. What exactly we compute
    // here is described in more detail in the documentation of that class.
    Vector<float> estimated_error_per_cell(
      this->triangulation->n_active_cells());
    std::map<types::boundary_id, const Function<dim> *> dummy_function_map;
    KellyErrorEstimator<dim>::estimate(this->dof_handler,
                                       *this->face_quadrature,
                                       dummy_function_map,
                                       this->solution,
                                       estimated_error_per_cell);

    // Next weigh each entry in the vector of indicators by the value of the
    // function given to the constructor, evaluated at the cell center. We
    // need to write the result into the vector entry that corresponds to the
    // current cell, which we can obtain by asking the cell what its index
    // among all active cells is using CellAccessor::active_cell_index(). (In
    // reality, this index is zero for the first cell we handle in the loop,
    // one for the second cell, etc., and we could as well just keep track of
    // this index using an integer counter; but using
    // CellAccessor::active_cell_index() makes this more explicit.)
    for (const auto &cell : this->dof_handler.active_cell_iterators())
      estimated_error_per_cell(cell->active_cell_index()) *=
        weighting_function->value(cell->center());

    GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);
    this->triangulation->execute_coarsening_and_refinement();
  }


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

  template <int dim>
  DualSolver<dim>::DualSolver(
    Triangulation<dim> &                           triangulation,
    const FiniteElement<dim> &                     fe,
    const Quadrature<dim> &                        quadrature,
    const Quadrature<dim - 1> &                    face_quadrature,
    const DualFunctional::DualFunctionalBase<dim> &dual_functional)
    : Base<dim>(triangulation)
    , Solver<dim>(triangulation,
                  fe,
                  quadrature,
                  face_quadrature,
                  boundary_values)
    , dual_functional(&dual_functional)
  {}



  template <int dim>
  void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const
  {
    dual_functional->assemble_rhs(this->dof_handler, rhs);
  }


  // @sect4{The WeightedResidual class}

  // Here finally comes the main class of this program, the one that
  // implements the dual weighted residual error estimator. It joins the
  // primal and dual solver classes to use them for the computation of
  // primal and dual solutions, and implements the error representation
  // formula for use as error estimate and mesh refinement.
  //
  // The first few of the functions of this class are mostly overriders of
  // the respective functions of the base class:
  template <int dim>
  class WeightedResidual : public PrimalSolver<dim>, public DualSolver<dim>
  {
  public:
    WeightedResidual(
      Triangulation<dim> &                           coarse_grid,
      const FiniteElement<dim> &                     primal_fe,
      const FiniteElement<dim> &                     dual_fe,
      const Quadrature<dim> &                        quadrature,
      const Quadrature<dim - 1> &                    face_quadrature,
      const Function<dim> &                          rhs_function,
      const Function<dim> &                          boundary_values,
      const DualFunctional::DualFunctionalBase<dim> &dual_functional);

    virtual void solve_problem() override;

    virtual void postprocess(
      const Evaluation::EvaluationBase<dim> &postprocessor) const override;

    virtual unsigned int n_dofs() const override;

    virtual void refine_grid() override;

    virtual void output_solution() const override;

  private:
    // In the private section, we have two functions that are used to call
    // the <code>solve_problem</code> functions of the primal and dual base
    // classes. These two functions will be called in parallel by the
    // <code>solve_problem</code> function of this class.
    void solve_primal_problem();
    void solve_dual_problem();
    // Then declare abbreviations for active cell iterators, to avoid that
    // we have to write this lengthy name over and over again:

    using active_cell_iterator =
      typename DoFHandler<dim>::active_cell_iterator;

    // Next, declare a data type that we will us to store the contribution
    // of faces to the error estimator. The idea is that we can compute the
    // face terms from each of the two cells to this face, as they are the
    // same when viewed from both sides. What we will do is to compute them
    // only once, based on some rules explained below which of the two
    // adjacent cells will be in charge to do so. We then store the
    // contribution of each face in a map mapping faces to their values, and
    // only collect the contributions for each cell by looping over the
    // cells a second time and grabbing the values from the map.
    //
    // The data type of this map is declared here:
    using FaceIntegrals =
      typename std::map<typename DoFHandler<dim>::face_iterator, double>;

    // In the computation of the error estimates on cells and faces, we need
    // a number of helper objects, such as <code>FEValues</code> and
    // <code>FEFaceValues</code> functions, but also temporary objects
    // storing the values and gradients of primal and dual solutions, for
    // example. These fields are needed in the three functions that do the
    // integration on cells, and regular and irregular faces, respectively.
    //
    // There are three reasonable ways to provide these fields: first, as
    // local variables in the function that needs them; second, as member
    // variables of this class; third, as arguments passed to that function.
    //
    // These three alternatives all have drawbacks: the third that their
    // number is not negligible and would make calling these functions a
    // lengthy enterprise. The second has the drawback that it disallows
    // parallelization, since the threads that will compute the error
    // estimate have to have their own copies of these variables each, so
    // member variables of the enclosing class will not work. The first
    // approach, although straightforward, has a subtle but important
    // drawback: we will call these functions over and over again, many
    // thousands of times maybe; it now turns out that allocating
    // vectors and other objects that need memory from the heap is an
    // expensive business in terms of run-time, since memory allocation is
    // expensive when several threads are involved. It is thus
    // significantly better to allocate the memory only once, and recycle
    // the objects as often as possible.
    //
    // What to do? Our answer is to use a variant of the third strategy.
    // In fact, this is exactly what the WorkStream concept is supposed to
    // do (we have already introduced it above, but see also @ref threads).
    // To avoid that we have to give these functions a dozen or so
    // arguments, we pack all these variables into two structures, one which
    // is used for the computations on cells, the other doing them on the
    // faces. Both are then joined into the WeightedResidualScratchData class
    // that will serve as the "scratch data" class of the WorkStream concept:
    struct CellData
    {
      FEValues<dim>                           fe_values;
      const SmartPointer<const Function<dim>> right_hand_side;

      std::vector<double> cell_residual;
      std::vector<double> rhs_values;
      std::vector<double> dual_weights;
      std::vector<double> cell_laplacians;
      CellData(const FiniteElement<dim> &fe,
               const Quadrature<dim> &   quadrature,
               const Function<dim> &     right_hand_side);
      CellData(const CellData &cell_data);
    };

    struct FaceData
    {
      FEFaceValues<dim>    fe_face_values_cell;
      FEFaceValues<dim>    fe_face_values_neighbor;
      FESubfaceValues<dim> fe_subface_values_cell;

      std::vector<double>                  jump_residual;
      std::vector<double>                  dual_weights;
      typename std::vector<Tensor<1, dim>> cell_grads;
      typename std::vector<Tensor<1, dim>> neighbor_grads;
      FaceData(const FiniteElement<dim> & fe,
               const Quadrature<dim - 1> &face_quadrature);
      FaceData(const FaceData &face_data);
    };

    struct WeightedResidualScratchData
    {
      WeightedResidualScratchData(
        const FiniteElement<dim> & primal_fe,
        const Quadrature<dim> &    primal_quadrature,
        const Quadrature<dim - 1> &primal_face_quadrature,
        const Function<dim> &      rhs_function,
        const Vector<double> &     primal_solution,
        const Vector<double> &     dual_weights);

      WeightedResidualScratchData(
        const WeightedResidualScratchData &scratch_data);

      CellData       cell_data;
      FaceData       face_data;
      Vector<double> primal_solution;
      Vector<double> dual_weights;
    };


    // WorkStream::run generally wants both a scratch object and a copy
    // object. Here, for reasons similar to what we had in step-9 when
    // discussing the computation of an approximation of the gradient, we
    // don't actually need a "copy data" structure. Since WorkStream insists
    // on having one of these, we just declare an empty structure that does
    // nothing other than being there.
    struct WeightedResidualCopyData
    {};



    // Regarding the evaluation of the error estimator, we have one driver
    // function that uses WorkStream::run() to call the second function on
    // every cell:
    void estimate_error(Vector<float> &error_indicators) const;

    void estimate_on_one_cell(const active_cell_iterator & cell,
                              WeightedResidualScratchData &scratch_data,
                              WeightedResidualCopyData &   copy_data,
                              Vector<float> &              error_indicators,
                              FaceIntegrals &face_integrals) const;

    // Then we have functions that do the actual integration of the error
    // representation formula. They will treat the terms on the cell
    // interiors, on those faces that have no hanging nodes, and on those
    // faces with hanging nodes, respectively:
    void integrate_over_cell(const active_cell_iterator &cell,
                             const Vector<double> &      primal_solution,
                             const Vector<double> &      dual_weights,
                             CellData &                  cell_data,
                             Vector<float> &error_indicators) const;

    void integrate_over_regular_face(const active_cell_iterator &cell,
                                     const unsigned int          face_no,
                                     const Vector<double> &primal_solution,
                                     const Vector<double> &dual_weights,
                                     FaceData &            face_data,
                                     FaceIntegrals &face_integrals) const;
    void integrate_over_irregular_face(const active_cell_iterator &cell,
                                       const unsigned int          face_no,
                                       const Vector<double> &primal_solution,
                                       const Vector<double> &dual_weights,
                                       FaceData &            face_data,
                                       FaceIntegrals &face_integrals) const;
  };



  // In the implementation of this class, we first have the constructors of
  // the <code>CellData</code> and <code>FaceData</code> member classes, and
  // the <code>WeightedResidual</code> constructor. They only initialize
  // fields to their correct lengths, so we do not have to discuss them in
  // too much detail:
  template <int dim>
  WeightedResidual<dim>::CellData::CellData(
    const FiniteElement<dim> &fe,
    const Quadrature<dim> &   quadrature,
    const Function<dim> &     right_hand_side)
    : fe_values(fe,
                quadrature,
                update_values | update_hessians | update_quadrature_points |
                  update_JxW_values)
    , right_hand_side(&right_hand_side)
    , cell_residual(quadrature.size())
    , rhs_values(quadrature.size())
    , dual_weights(quadrature.size())
    , cell_laplacians(quadrature.size())
  {}



  template <int dim>
  WeightedResidual<dim>::CellData::CellData(const CellData &cell_data)
    : fe_values(cell_data.fe_values.get_fe(),
                cell_data.fe_values.get_quadrature(),
                update_values | update_hessians | update_quadrature_points |
                  update_JxW_values)
    , right_hand_side(cell_data.right_hand_side)
    , cell_residual(cell_data.cell_residual)
    , rhs_values(cell_data.rhs_values)
    , dual_weights(cell_data.dual_weights)
    , cell_laplacians(cell_data.cell_laplacians)
  {}



  template <int dim>
  WeightedResidual<dim>::FaceData::FaceData(
    const FiniteElement<dim> & fe,
    const Quadrature<dim - 1> &face_quadrature)
    : fe_face_values_cell(fe,
                          face_quadrature,
                          update_values | update_gradients |
                            update_JxW_values | update_normal_vectors)
    , fe_face_values_neighbor(fe,
                              face_quadrature,
                              update_values | update_gradients |
                                update_JxW_values | update_normal_vectors)
    , fe_subface_values_cell(fe, face_quadrature, update_gradients)
  {
    const unsigned int n_face_q_points = face_quadrature.size();

    jump_residual.resize(n_face_q_points);
    dual_weights.resize(n_face_q_points);
    cell_grads.resize(n_face_q_points);
    neighbor_grads.resize(n_face_q_points);
  }



  template <int dim>
  WeightedResidual<dim>::FaceData::FaceData(const FaceData &face_data)
    : fe_face_values_cell(face_data.fe_face_values_cell.get_fe(),
                          face_data.fe_face_values_cell.get_quadrature(),
                          update_values | update_gradients |
                            update_JxW_values | update_normal_vectors)
    , fe_face_values_neighbor(
        face_data.fe_face_values_neighbor.get_fe(),
        face_data.fe_face_values_neighbor.get_quadrature(),
        update_values | update_gradients | update_JxW_values |
          update_normal_vectors)
    , fe_subface_values_cell(
        face_data.fe_subface_values_cell.get_fe(),
        face_data.fe_subface_values_cell.get_quadrature(),
        update_gradients)
    , jump_residual(face_data.jump_residual)
    , dual_weights(face_data.dual_weights)
    , cell_grads(face_data.cell_grads)
    , neighbor_grads(face_data.neighbor_grads)
  {}



  template <int dim>
  WeightedResidual<dim>::WeightedResidualScratchData::
    WeightedResidualScratchData(
      const FiniteElement<dim> & primal_fe,
      const Quadrature<dim> &    primal_quadrature,
      const Quadrature<dim - 1> &primal_face_quadrature,
      const Function<dim> &      rhs_function,
      const Vector<double> &     primal_solution,
      const Vector<double> &     dual_weights)
    : cell_data(primal_fe, primal_quadrature, rhs_function)
    , face_data(primal_fe, primal_face_quadrature)
    , primal_solution(primal_solution)
    , dual_weights(dual_weights)
  {}

  template <int dim>
  WeightedResidual<dim>::WeightedResidualScratchData::
    WeightedResidualScratchData(
      const WeightedResidualScratchData &scratch_data)
    : cell_data(scratch_data.cell_data)
    , face_data(scratch_data.face_data)
    , primal_solution(scratch_data.primal_solution)
    , dual_weights(scratch_data.dual_weights)
  {}



  template <int dim>
  WeightedResidual<dim>::WeightedResidual(
    Triangulation<dim> &                           coarse_grid,
    const FiniteElement<dim> &                     primal_fe,
    const FiniteElement<dim> &                     dual_fe,
    const Quadrature<dim> &                        quadrature,
    const Quadrature<dim - 1> &                    face_quadrature,
    const Function<dim> &                          rhs_function,
    const Function<dim> &                          bv,
    const DualFunctional::DualFunctionalBase<dim> &dual_functional)
    : Base<dim>(coarse_grid)
    , PrimalSolver<dim>(coarse_grid,
                        primal_fe,
                        quadrature,
                        face_quadrature,
                        rhs_function,
                        bv)
    , DualSolver<dim>(coarse_grid,
                      dual_fe,
                      quadrature,
                      face_quadrature,
                      dual_functional)
  {}


  // The next five functions are boring, as they simply relay their work to
  // the base classes. The first calls the primal and dual solvers in
  // parallel, while postprocessing the solution and retrieving the number
  // of degrees of freedom is done by the primal class.
  template <int dim>
  void WeightedResidual<dim>::solve_problem()
  {
    Threads::TaskGroup<void> tasks;
    tasks +=
      Threads::new_task(&WeightedResidual<dim>::solve_primal_problem, *this);
    tasks +=
      Threads::new_task(&WeightedResidual<dim>::solve_dual_problem, *this);
    tasks.join_all();
  }


  template <int dim>
  void WeightedResidual<dim>::solve_primal_problem()
  {
    PrimalSolver<dim>::solve_problem();
  }

  template <int dim>
  void WeightedResidual<dim>::solve_dual_problem()
  {
    DualSolver<dim>::solve_problem();
  }


  template <int dim>
  void WeightedResidual<dim>::postprocess(
    const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    PrimalSolver<dim>::postprocess(postprocessor);
  }


  template <int dim>
  unsigned int WeightedResidual<dim>::n_dofs() const
  {
    return PrimalSolver<dim>::n_dofs();
  }



  // Now, it is becoming more interesting: the <code>refine_grid()</code>
  // function asks the error estimator to compute the cell-wise error
  // indicators, then uses their absolute values for mesh refinement.
  template <int dim>
  void WeightedResidual<dim>::refine_grid()
  {
    // First call the function that computes the cell-wise and global error:
    Vector<float> error_indicators(this->triangulation->n_active_cells());
    estimate_error(error_indicators);

    // Then note that marking cells for refinement or coarsening only works
    // if all indicators are positive, to allow their comparison. Thus, drop
    // the signs on all these indicators:
    for (float &error_indicator : error_indicators)
      error_indicator = std::fabs(error_indicator);

    // Finally, we can select between different strategies for
    // refinement. The default here is to refine those cells with the
    // largest error indicators that make up for a total of 80 per cent of
    // the error, while we coarsen those with the smallest indicators that
    // make up for the bottom 2 per cent of the error.
    GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                      error_indicators,
                                                      0.8,
                                                      0.02);
    this->triangulation->execute_coarsening_and_refinement();
  }


  // Since we want to output both the primal and the dual solution, we
  // overload the <code>output_solution</code> function. The only
  // interesting feature of this function is that the primal and dual
  // solutions are defined on different finite element spaces, which is not
  // the format the <code>DataOut</code> class expects. Thus, we have to
  // transfer them to a common finite element space. Since we want the
  // solutions only to see them qualitatively, we contend ourselves with
  // interpolating the dual solution to the (smaller) primal space. For the
  // interpolation, there is a library function, that takes a
  // AffineConstraints object including the hanging node
  // constraints. The rest is standard.
  template <int dim>
  void WeightedResidual<dim>::output_solution() const
  {
    AffineConstraints<double> primal_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler,
                                            primal_hanging_node_constraints);
    primal_hanging_node_constraints.close();
    Vector<double> dual_solution(PrimalSolver<dim>::dof_handler.n_dofs());
    FETools::interpolate(DualSolver<dim>::dof_handler,
                         DualSolver<dim>::solution,
                         PrimalSolver<dim>::dof_handler,
                         primal_hanging_node_constraints,
                         dual_solution);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(PrimalSolver<dim>::dof_handler);

    // Add the data vectors for which we want output. Add them both, the
    // <code>DataOut</code> functions can handle as many data vectors as you
    // wish to write to output:
    data_out.add_data_vector(PrimalSolver<dim>::solution, "primal_solution");
    data_out.add_data_vector(dual_solution, "dual_solution");

    data_out.build_patches();

    std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                      ".vtu");
    data_out.write(out, DataOutBase::vtu);
  }


  // @sect3{Estimating errors}

  // @sect4{Error estimation driver functions}
  //
  // As for the actual computation of error estimates, let's start with the
  // function that drives all this, i.e. calls those functions that actually
  // do the work, and finally collects the results.
  template <int dim>
  void
  WeightedResidual<dim>::estimate_error(Vector<float> &error_indicators) const
  {
    // The first task in computing the error is to set up vectors that
    // denote the primal solution, and the weights (z-z_h)=(z-I_hz), both in
    // the finite element space for which we have computed the dual
    // solution. For this, we have to interpolate the primal solution to the
    // dual finite element space, and to subtract the interpolation of the
    // computed dual solution to the primal finite element
    // space. Fortunately, the library provides functions for the
    // interpolation into larger or smaller finite element spaces, so this
    // is mostly obvious.
    //
    // First, let's do that for the primal solution: it is cell-wise
    // interpolated into the finite element space in which we have solved
    // the dual problem: But, again as in the
    // <code>WeightedResidual::output_solution</code> function we first need
    // to create an AffineConstraints object including the hanging node
    // constraints, but this time of the dual finite element space.
    AffineConstraints<double> dual_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler,
                                            dual_hanging_node_constraints);
    dual_hanging_node_constraints.close();
    Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs());
    FETools::interpolate(PrimalSolver<dim>::dof_handler,
                         PrimalSolver<dim>::solution,
                         DualSolver<dim>::dof_handler,
                         dual_hanging_node_constraints,
                         primal_solution);

    // Then for computing the interpolation of the numerically approximated
    // dual solution z into the finite element space of the primal solution
    // and subtracting it from z: use the
    // <code>interpolate_difference</code> function, that gives (z-I_hz) in
    // the element space of the dual solution.
    AffineConstraints<double> primal_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler,
                                            primal_hanging_node_constraints);
    primal_hanging_node_constraints.close();
    Vector<double> dual_weights(DualSolver<dim>::dof_handler.n_dofs());
    FETools::interpolation_difference(DualSolver<dim>::dof_handler,
                                      dual_hanging_node_constraints,
                                      DualSolver<dim>::solution,
                                      PrimalSolver<dim>::dof_handler,
                                      primal_hanging_node_constraints,
                                      dual_weights);

    // Note that this could probably have been more efficient since those
    // constraints have been used previously when assembling matrix and
    // right hand side for the primal problem and writing out the dual
    // solution. We leave the optimization of the program in this respect as
    // an exercise.

    // Having computed the dual weights we now proceed with computing the
    // cell and face residuals of the primal solution. First we set up a map
    // between face iterators and their jump term contributions of faces to
    // the error estimator. The reason is that we compute the jump terms
    // only once, from one side of the face, and want to collect them only
    // afterwards when looping over all cells a second time.
    //
    // We initialize this map already with a value of -1e20 for all faces,
    // since this value will stand out in the results if something should go
    // wrong and we fail to compute the value for a face for some
    // reason. Secondly, this initialization already makes the std::map
    // object allocate all objects it may possibly need. This is important
    // since we will write into this structure from parallel threads,
    // and doing so would not be thread-safe if the map needed to allocate
    // memory and thereby reshape its data structures. In other words, the
    // initial initialization relieves us from the necessity to synchronize
    // the threads through a mutex each time they write to (and modify the
    // structure of) this map.
    FaceIntegrals face_integrals;
    for (const auto &cell :
         DualSolver<dim>::dof_handler.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
        face_integrals[face] = -1e20;

    auto worker = [this,
                   &error_indicators,
                   &face_integrals](const active_cell_iterator & cell,
                                    WeightedResidualScratchData &scratch_data,
                                    WeightedResidualCopyData &   copy_data) {
      this->estimate_on_one_cell(
        cell, scratch_data, copy_data, error_indicators, face_integrals);
    };

    auto do_nothing_copier =
      std::function<void(const WeightedResidualCopyData &)>();

    // Then hand it all off to WorkStream::run() to compute the
    // estimators for all cells in parallel:
    WorkStream::run(
      DualSolver<dim>::dof_handler.begin_active(),
      DualSolver<dim>::dof_handler.end(),
      worker,
      do_nothing_copier,
      WeightedResidualScratchData(*DualSolver<dim>::fe,
                                  *DualSolver<dim>::quadrature,
                                  *DualSolver<dim>::face_quadrature,
                                  *this->rhs_function,
                                  primal_solution,
                                  dual_weights),
      WeightedResidualCopyData());

    // Once the error contributions are computed, sum them up. For this,
    // note that the cell terms are already set, and that only the edge
    // terms need to be collected. Thus, loop over all cells and their
    // faces, make sure that the contributions of each of the faces are
    // there, and add them up. Only take minus one half of the jump term,
    // since the other half will be taken by the neighboring cell.
    unsigned int present_cell = 0;
    for (const auto &cell :
         DualSolver<dim>::dof_handler.active_cell_iterators())
      {
        for (const auto &face : cell->face_iterators())
          {
            Assert(face_integrals.find(face) != face_integrals.end(),
                   ExcInternalError());
            error_indicators(present_cell) -= 0.5 * face_integrals[face];
          }
        ++present_cell;
      }
    std::cout << "   Estimated error="
              << std::accumulate(error_indicators.begin(),
                                 error_indicators.end(),
                                 0.)
              << std::endl;
  }


  // @sect4{Estimating on a single cell}

  // Next we have the function that is called to estimate the error on a
  // single cell. The function may be called multiple times if the library was
  // configured to use multithreading. Here it goes:
  template <int dim>
  void WeightedResidual<dim>::estimate_on_one_cell(
    const active_cell_iterator & cell,
    WeightedResidualScratchData &scratch_data,
    WeightedResidualCopyData &   copy_data,
    Vector<float> &              error_indicators,
    FaceIntegrals &              face_integrals) const
  {
    // Because of WorkStream, estimate_on_one_cell requires a CopyData object
    // even if it is no used. The next line silences a warning about this
    // unused variable.
    (void)copy_data;

    // First task on each cell is to compute the cell residual
    // contributions of this cell, and put them into the
    // <code>error_indicators</code> variable:
    integrate_over_cell(cell,
                        scratch_data.primal_solution,
                        scratch_data.dual_weights,
                        scratch_data.cell_data,
                        error_indicators);

    // After computing the cell terms, turn to the face terms. For this,
    // loop over all faces of the present cell, and see whether
    // something needs to be computed on it:
    for (const auto face_no : cell->face_indices())
      {
        // First, if this face is part of the boundary, then there is
        // nothing to do. However, to make things easier when summing up
        // the contributions of the faces of cells, we enter this face
        // into the list of faces with a zero contribution to the error.
        if (cell->face(face_no)->at_boundary())
          {
            face_integrals[cell->face(face_no)] = 0;
            continue;
          }

        // Next, note that since we want to compute the jump terms on
        // each face only once although we access it twice (if it is not
        // at the boundary), we have to define some rules who is
        // responsible for computing on a face:
        //
        // First, if the neighboring cell is on the same level as this
        // one, i.e. neither further refined not coarser, then the one
        // with the lower index within this level does the work. In
        // other words: if the other one has a lower index, then skip
        // work on this face:
        if ((cell->neighbor(face_no)->has_children() == false) &&
            (cell->neighbor(face_no)->level() == cell->level()) &&
            (cell->neighbor(face_no)->index() < cell->index()))
          continue;

        // Likewise, we always work from the coarser cell if this and
        // its neighbor differ in refinement. Thus, if the neighboring
        // cell is less refined than the present one, then do nothing
        // since we integrate over the subfaces when we visit the coarse
        // cell.
        if (cell->at_boundary(face_no) == false)
          if (cell->neighbor(face_no)->level() < cell->level())
            continue;


        // Now we know that we are in charge here, so actually compute
        // the face jump terms. If the face is a regular one, i.e.  the
        // other side's cell is neither coarser not finer than this
        // cell, then call one function, and if the cell on the other
        // side is further refined, then use another function. Note that
        // the case that the cell on the other side is coarser cannot
        // happen since we have decided above that we handle this case
        // when we pass over that other cell.
        if (cell->face(face_no)->has_children() == false)
          integrate_over_regular_face(cell,
                                      face_no,
                                      scratch_data.primal_solution,
                                      scratch_data.dual_weights,
                                      scratch_data.face_data,
                                      face_integrals);
        else
          integrate_over_irregular_face(cell,
                                        face_no,
                                        scratch_data.primal_solution,
                                        scratch_data.dual_weights,
                                        scratch_data.face_data,
                                        face_integrals);
      }
  }


  // @sect4{Computing cell term error contributions}

  // As for the actual computation of the error contributions, first turn to
  // the cell terms:
  template <int dim>
  void WeightedResidual<dim>::integrate_over_cell(
    const active_cell_iterator &cell,
    const Vector<double> &      primal_solution,
    const Vector<double> &      dual_weights,
    CellData &                  cell_data,
    Vector<float> &             error_indicators) const
  {
    // The tasks to be done are what appears natural from looking at the
    // error estimation formula: first get the right hand side and Laplacian
    // of the numerical solution at the quadrature points for the cell
    // residual,
    cell_data.fe_values.reinit(cell);
    cell_data.right_hand_side->value_list(
      cell_data.fe_values.get_quadrature_points(), cell_data.rhs_values);
    cell_data.fe_values.get_function_laplacians(primal_solution,
                                                cell_data.cell_laplacians);

    // ...then get the dual weights...
    cell_data.fe_values.get_function_values(dual_weights,
                                            cell_data.dual_weights);

    // ...and finally build the sum over all quadrature points and store it
    // with the present cell:
    double sum = 0;
    for (unsigned int p = 0; p < cell_data.fe_values.n_quadrature_points; ++p)
      sum += ((cell_data.rhs_values[p] + cell_data.cell_laplacians[p]) *
              cell_data.dual_weights[p] * cell_data.fe_values.JxW(p));
    error_indicators(cell->active_cell_index()) += sum;
  }


  // @sect4{Computing edge term error contributions -- 1}

  // On the other hand, computation of the edge terms for the error estimate
  // is not so simple. First, we have to distinguish between faces with and
  // without hanging nodes. Because it is the simple case, we first consider
  // the case without hanging nodes on a face (let's call this the `regular'
  // case):
  template <int dim>
  void WeightedResidual<dim>::integrate_over_regular_face(
    const active_cell_iterator &cell,
    const unsigned int          face_no,
    const Vector<double> &      primal_solution,
    const Vector<double> &      dual_weights,
    FaceData &                  face_data,
    FaceIntegrals &             face_integrals) const
  {
    const unsigned int n_q_points =
      face_data.fe_face_values_cell.n_quadrature_points;

    // The first step is to get the values of the gradients at the
    // quadrature points of the finite element field on the present
    // cell. For this, initialize the <code>FEFaceValues</code> object
    // corresponding to this side of the face, and extract the gradients
    // using that object.
    face_data.fe_face_values_cell.reinit(cell, face_no);
    face_data.fe_face_values_cell.get_function_gradients(
      primal_solution, face_data.cell_grads);

    // The second step is then to extract the gradients of the finite
    // element solution at the quadrature points on the other side of the
    // face, i.e. from the neighboring cell.
    //
    // For this, do a sanity check before: make sure that the neighbor
    // actually exists (yes, we should not have come here if the neighbor
    // did not exist, but in complicated software there are bugs, so better
    // check this), and if this is not the case throw an error.
    Assert(cell->neighbor(face_no).state() == IteratorState::valid,
           ExcInternalError());
    // If we have that, then we need to find out with which face of the
    // neighboring cell we have to work, i.e. the <code>how-many'th</code> the
    // neighbor the present cell is of the cell behind the present face. For
    // this, there is a function, and we put the result into a variable with
    // the name <code>neighbor_neighbor</code>:
    const unsigned int neighbor_neighbor =
      cell->neighbor_of_neighbor(face_no);
    // Then define an abbreviation for the neighbor cell, initialize the
    // <code>FEFaceValues</code> object on that cell, and extract the
    // gradients on that cell:
    const active_cell_iterator neighbor = cell->neighbor(face_no);
    face_data.fe_face_values_neighbor.reinit(neighbor, neighbor_neighbor);
    face_data.fe_face_values_neighbor.get_function_gradients(
      primal_solution, face_data.neighbor_grads);

    // Now that we have the gradients on this and the neighboring cell,
    // compute the jump residual by multiplying the jump in the gradient
    // with the normal vector:
    for (unsigned int p = 0; p < n_q_points; ++p)
      face_data.jump_residual[p] =
        ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) *
         face_data.fe_face_values_cell.normal_vector(p));

    // Next get the dual weights for this face:
    face_data.fe_face_values_cell.get_function_values(dual_weights,
                                                      face_data.dual_weights);

    // Finally, we have to compute the sum over jump residuals, dual
    // weights, and quadrature weights, to get the result for this face:
    double face_integral = 0;
    for (unsigned int p = 0; p < n_q_points; ++p)
      face_integral +=
        (face_data.jump_residual[p] * face_data.dual_weights[p] *
         face_data.fe_face_values_cell.JxW(p));

    // Double check that the element already exists and that it was not
    // already written to...
    Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(),
           ExcInternalError());
    Assert(face_integrals[cell->face(face_no)] == -1e20, ExcInternalError());

    // ...then store computed value at assigned location. Note that the
    // stored value does not contain the factor 1/2 that appears in the
    // error representation. The reason is that the term actually does not
    // have this factor if we loop over all faces in the triangulation, but
    // only appears if we write it as a sum over all cells and all faces of
    // each cell; we thus visit the same face twice. We take account of this
    // by using this factor -1/2 later, when we sum up the contributions for
    // each cell individually.
    face_integrals[cell->face(face_no)] = face_integral;
  }


  // @sect4{Computing edge term error contributions -- 2}

  // We are still missing the case of faces with hanging nodes. This is what
  // is covered in this function:
  template <int dim>
  void WeightedResidual<dim>::integrate_over_irregular_face(
    const active_cell_iterator &cell,
    const unsigned int          face_no,
    const Vector<double> &      primal_solution,
    const Vector<double> &      dual_weights,
    FaceData &                  face_data,
    FaceIntegrals &             face_integrals) const
  {
    // First again two abbreviations, and some consistency checks whether
    // the function is called only on faces for which it is supposed to be
    // called:
    const unsigned int n_q_points =
      face_data.fe_face_values_cell.n_quadrature_points;

    const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
    const typename DoFHandler<dim>::cell_iterator neighbor =
      cell->neighbor(face_no);
    Assert(neighbor.state() == IteratorState::valid, ExcInternalError());
    Assert(neighbor->has_children(), ExcInternalError());
    (void)neighbor;

    // Then find out which neighbor the present cell is of the adjacent
    // cell. Note that we will operate on the children of this adjacent
    // cell, but that their orientation is the same as that of their mother,
    // i.e. the neighbor direction is the same.
    const unsigned int neighbor_neighbor =
      cell->neighbor_of_neighbor(face_no);

    // Then simply do everything we did in the previous function for one
    // face for all the sub-faces now:
    for (unsigned int subface_no = 0; subface_no < face->n_children();
         ++subface_no)
      {
        // Start with some checks again: get an iterator pointing to the
        // cell behind the present subface and check whether its face is a
        // subface of the one we are considering. If that were not the case,
        // then there would be either a bug in the
        // <code>neighbor_neighbor</code> function called above, or -- worse
        // -- some function in the library did not keep to some underlying
        // assumptions about cells, their children, and their faces. In any
        // case, even though this assertion should not be triggered, it does
        // not harm to be cautious, and in optimized mode computations the
        // assertion will be removed anyway.
        const active_cell_iterator neighbor_child =
          cell->neighbor_child_on_subface(face_no, subface_no);
        Assert(neighbor_child->face(neighbor_neighbor) ==
                 cell->face(face_no)->child(subface_no),
               ExcInternalError());

        // Now start the work by again getting the gradient of the solution
        // first at this side of the interface,
        face_data.fe_subface_values_cell.reinit(cell, face_no, subface_no);
        face_data.fe_subface_values_cell.get_function_gradients(
          primal_solution, face_data.cell_grads);
        // then at the other side,
        face_data.fe_face_values_neighbor.reinit(neighbor_child,
                                                 neighbor_neighbor);
        face_data.fe_face_values_neighbor.get_function_gradients(
          primal_solution, face_data.neighbor_grads);

        // and finally building the jump residuals. Since we take the normal
        // vector from the other cell this time, revert the sign of the
        // first term compared to the other function:
        for (unsigned int p = 0; p < n_q_points; ++p)
          face_data.jump_residual[p] =
            ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) *
             face_data.fe_face_values_neighbor.normal_vector(p));

        // Then get dual weights:
        face_data.fe_face_values_neighbor.get_function_values(
          dual_weights, face_data.dual_weights);

        // At last, sum up the contribution of this sub-face, and set it in
        // the global map:
        double face_integral = 0;
        for (unsigned int p = 0; p < n_q_points; ++p)
          face_integral +=
            (face_data.jump_residual[p] * face_data.dual_weights[p] *
             face_data.fe_face_values_neighbor.JxW(p));
        face_integrals[neighbor_child->face(neighbor_neighbor)] =
          face_integral;
      }

    // Once the contributions of all sub-faces are computed, loop over all
    // sub-faces to collect and store them with the mother face for simple
    // use when later collecting the error terms of cells. Again make safety
    // checks that the entries for the sub-faces have been computed and do
    // not carry an invalid value.
    double sum = 0;
    for (unsigned int subface_no = 0; subface_no < face->n_children();
         ++subface_no)
      {
        Assert(face_integrals.find(face->child(subface_no)) !=
                 face_integrals.end(),
               ExcInternalError());
        Assert(face_integrals[face->child(subface_no)] != -1e20,
               ExcInternalError());

        sum += face_integrals[face->child(subface_no)];
      }
    // Finally store the value with the parent face.
    face_integrals[face] = sum;
  }




  } // namespace LaplaceSolver
}

#endif