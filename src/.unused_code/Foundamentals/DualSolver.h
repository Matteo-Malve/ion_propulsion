#ifndef ION_PROPULSION_DUALSOLVER_H
#define ION_PROPULSION_DUALSOLVER_H

#include "DualFunctional.h"
#include "Solver.h"
#include "HelperFunctions.h"

template <int dim>
class DualSolver : public Solver<dim>
{
public:
    // Constructor
    DualSolver(
            Triangulation<dim> &                           triangulation,
            const FiniteElement<dim> &                     fe,
            const Quadrature<dim> &                        quadrature,
            const Quadrature<dim - 1> &                    face_quadrature,
            const DualFunctionalBase<dim> &dual_functional);
    // Override solve method
    virtual void solve_problem() override;
    virtual void output_solution() override;
protected:
    // DualFunctional class comes into play
    const SmartPointer<const DualFunctionalBase<dim>>
            dual_functional;
    virtual void assemble_rhs(Vector<double> &rhs) const override;
    virtual void apply_boundary_conditions() override;
};

// CONSTRUCTOR
template <int dim>
DualSolver<dim>::DualSolver(
        Triangulation<dim> &                           triangulation_,
        const FiniteElement<dim> &                     fe_,
        const Quadrature<dim> &                        quadrature_,
        const Quadrature<dim - 1> &                    face_quadrature_,
        const DualFunctionalBase<dim> &dual_functional_)
        : Base<dim>(triangulation_)
        , Solver<dim>(triangulation_,
                      fe_,
                      quadrature_,
                      face_quadrature_)
        , dual_functional(&dual_functional_)
{}

// ASSEMBLE_rhs override
template <int dim>
void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const
{
    dual_functional->assemble_rhs(this->dof_handler, rhs);

}

template<int dim>
void DualSolver<dim>::apply_boundary_conditions() {
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    // Associate BCs to corresponding Boundary Ids. The latter are directly set in Gmsh, as physical curves
    VectorTools::interpolate_boundary_values(this->dof_handler,
                                             1, // Boundary id corresponding to Emitter
                                             Functions::ConstantFunction<dim>(0),
                                             emitter_boundary_values);

    VectorTools::interpolate_boundary_values(this->dof_handler,
                                             2,  // Boundary id corresponding to the Collector
                                             Functions::ConstantFunction<dim>(0),
                                             collector_boundary_values);
    // NOTE: z belongs to H1_O_GammaD
    //       By default, deal.ii will then set all non-specified boundary points to have homogeneous Neumann BCs.
    //       That's exactly what we need.

    MatrixTools::apply_boundary_values(emitter_boundary_values,
                                       this->system_matrix,
                                       this->solution,
                                       this->system_rhs);

    MatrixTools::apply_boundary_values(collector_boundary_values,
                                       this->system_matrix,
                                       this->solution,
                                       this->system_rhs);
}

template <int dim>
void DualSolver<dim>::solve_problem()
{
    this->setup_system();
    this->assemble_system();
    this->apply_boundary_conditions();
    this->solve_system();
}

template <int dim>
void DualSolver<dim>::output_solution()
{

    // Write solution to file. Format: .vtu
    GradientPostprocessor<dim> gradient_postprocessor;
    DataOut <dim> data_out;
    data_out.attach_dof_handler(this->dof_handler);
    data_out.add_data_vector(this->solution, "Potential");
    data_out.add_data_vector(this->solution, gradient_postprocessor);
    data_out.build_patches();
    std::ofstream output("dual_solution-" + std::to_string(this->refinement_cycle) + ".vtu");
    data_out.write_vtu(output);

}

#endif //ION_PROPULSION_DUALSOLVER_H
