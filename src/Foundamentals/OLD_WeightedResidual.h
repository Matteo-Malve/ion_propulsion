#ifndef GETPOT_WEIGHTEDRESIDUAL_H
#define GETPOT_WEIGHTEDRESIDUAL_H

#include "PrimalSolver.h"
#include "DualSolver.h"

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
            const DualFunctionalBase<dim> &dual_functional);

    virtual void solve_problem() override;

    virtual void postprocess(
            const Evaluation::EvaluationBase<dim> &postprocessor) const override;

    virtual unsigned int n_dofs() const override;

    virtual void refine_grid() override;

    virtual void output_solution() const override;

private:
    void solve_primal_problem();
    void solve_dual_problem();

    using active_cell_iterator =
            typename DoFHandler<dim>::active_cell_iterator;

    using FaceIntegrals =
            typename std::map<typename DoFHandler<dim>::face_iterator, double>;

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


    struct WeightedResidualCopyData
    {};



    void estimate_error(Vector<float> &error_indicators) const;

    void estimate_on_one_cell(const active_cell_iterator & cell,
                              WeightedResidualScratchData &scratch_data,
                              WeightedResidualCopyData &   copy_data,
                              Vector<float> &              error_indicators,
                              FaceIntegrals &face_integrals) const;

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



template <int dim>
void WeightedResidual<dim>::refine_grid()
{
    Vector<float> error_indicators(this->triangulation->n_active_cells());
    estimate_error(error_indicators);

    for (float &error_indicator : error_indicators)
        error_indicator = std::fabs(error_indicator);

    GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
                                                      error_indicators,
                                                      0.8,
                                                      0.02);
    this->triangulation->execute_coarsening_and_refinement();
}


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

    data_out.add_data_vector(PrimalSolver<dim>::solution, "primal_solution");
    data_out.add_data_vector(dual_solution, "dual_solution");

    data_out.build_patches();

    std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                      ".vtu");
    data_out.write(out, DataOutBase::vtu);
}



template <int dim>
void
WeightedResidual<dim>::estimate_error(Vector<float> &error_indicators) const
{
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



template <int dim>
void WeightedResidual<dim>::estimate_on_one_cell(
        const active_cell_iterator & cell,
        WeightedResidualScratchData &scratch_data,
        WeightedResidualCopyData &   copy_data,
        Vector<float> &              error_indicators,
        FaceIntegrals &              face_integrals) const
{
    (void)copy_data;

    integrate_over_cell(cell,
                        scratch_data.primal_solution,
                        scratch_data.dual_weights,
                        scratch_data.cell_data,
                        error_indicators);

    for (const auto face_no : cell->face_indices())
    {
        if (cell->face(face_no)->at_boundary())
        {
            face_integrals[cell->face(face_no)] = 0;
            continue;
        }

        if ((cell->neighbor(face_no)->has_children() == false) &&
            (cell->neighbor(face_no)->level() == cell->level()) &&
            (cell->neighbor(face_no)->index() < cell->index()))
            continue;

        if (cell->at_boundary(face_no) == false)
            if (cell->neighbor(face_no)->level() < cell->level())
                continue;


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



template <int dim>
void WeightedResidual<dim>::integrate_over_cell(
        const active_cell_iterator &cell,
        const Vector<double> &      primal_solution,
        const Vector<double> &      dual_weights,
        CellData &                  cell_data,
        Vector<float> &             error_indicators) const
{
    cell_data.fe_values.reinit(cell);
    cell_data.right_hand_side->value_list(
            cell_data.fe_values.get_quadrature_points(), cell_data.rhs_values);
    cell_data.fe_values.get_function_laplacians(primal_solution,
                                                cell_data.cell_laplacians);

    cell_data.fe_values.get_function_values(dual_weights,
                                            cell_data.dual_weights);

    double sum = 0;
    for (unsigned int p = 0; p < cell_data.fe_values.n_quadrature_points; ++p)
        sum += ((cell_data.rhs_values[p] + cell_data.cell_laplacians[p]) *
                cell_data.dual_weights[p] * cell_data.fe_values.JxW(p));
    error_indicators(cell->active_cell_index()) += sum;
}



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

    face_data.fe_face_values_cell.reinit(cell, face_no);
    face_data.fe_face_values_cell.get_function_gradients(
            primal_solution, face_data.cell_grads);

    Assert(cell->neighbor(face_no).state() == IteratorState::valid,
           ExcInternalError());
    const unsigned int neighbor_neighbor =
            cell->neighbor_of_neighbor(face_no);
    const active_cell_iterator neighbor = cell->neighbor(face_no);
    face_data.fe_face_values_neighbor.reinit(neighbor, neighbor_neighbor);
    face_data.fe_face_values_neighbor.get_function_gradients(
            primal_solution, face_data.neighbor_grads);

    for (unsigned int p = 0; p < n_q_points; ++p)
        face_data.jump_residual[p] =
                ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) *
                 face_data.fe_face_values_cell.normal_vector(p));

    face_data.fe_face_values_cell.get_function_values(dual_weights,
                                                      face_data.dual_weights);

    double face_integral = 0;
    for (unsigned int p = 0; p < n_q_points; ++p)
        face_integral +=
                (face_data.jump_residual[p] * face_data.dual_weights[p] *
                 face_data.fe_face_values_cell.JxW(p));

    Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(),
           ExcInternalError());
    Assert(face_integrals[cell->face(face_no)] == -1e20, ExcInternalError());

    face_integrals[cell->face(face_no)] = face_integral;
}



template <int dim>
void WeightedResidual<dim>::integrate_over_irregular_face(
        const active_cell_iterator &cell,
        const unsigned int          face_no,
        const Vector<double> &      primal_solution,
        const Vector<double> &      dual_weights,
        FaceData &                  face_data,
        FaceIntegrals &             face_integrals) const
{
    const unsigned int n_q_points =
            face_data.fe_face_values_cell.n_quadrature_points;

    const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
    const typename DoFHandler<dim>::cell_iterator neighbor =
            cell->neighbor(face_no);
    Assert(neighbor.state() == IteratorState::valid, ExcInternalError());
    Assert(neighbor->has_children(), ExcInternalError());
    (void)neighbor;

    const unsigned int neighbor_neighbor =
            cell->neighbor_of_neighbor(face_no);

    for (unsigned int subface_no = 0; subface_no < face->n_children();
         ++subface_no)
    {
        const active_cell_iterator neighbor_child =
                cell->neighbor_child_on_subface(face_no, subface_no);
        Assert(neighbor_child->face(neighbor_neighbor) ==
               cell->face(face_no)->child(subface_no),
               ExcInternalError());

        face_data.fe_subface_values_cell.reinit(cell, face_no, subface_no);
        face_data.fe_subface_values_cell.get_function_gradients(
                primal_solution, face_data.cell_grads);
        face_data.fe_face_values_neighbor.reinit(neighbor_child,
                                                 neighbor_neighbor);
        face_data.fe_face_values_neighbor.get_function_gradients(
                primal_solution, face_data.neighbor_grads);

        for (unsigned int p = 0; p < n_q_points; ++p)
            face_data.jump_residual[p] =
                    ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) *
                     face_data.fe_face_values_neighbor.normal_vector(p));

        face_data.fe_face_values_neighbor.get_function_values(
                dual_weights, face_data.dual_weights);

        double face_integral = 0;
        for (unsigned int p = 0; p < n_q_points; ++p)
            face_integral +=
                    (face_data.jump_residual[p] * face_data.dual_weights[p] *
                     face_data.fe_face_values_neighbor.JxW(p));
        face_integrals[neighbor_child->face(neighbor_neighbor)] =
                face_integral;
    }

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
    face_integrals[face] = sum;
}


#endif //GETPOT_WEIGHTEDRESIDUAL_H
