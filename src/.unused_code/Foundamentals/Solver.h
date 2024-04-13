#ifndef ION_PROPULSION_SOLVER_H
#define ION_PROPULSION_SOLVER_H
#include "Base.h"

template <int dim>
class Solver : public virtual Base<dim>
{
public:
    Solver(Triangulation<dim> &       triangulation,
           const FiniteElement<dim> & fe,
           const Quadrature<dim> &    quadrature,
           const Quadrature<dim - 1> &face_quadrature);

    virtual ~Solver() override;

    virtual unsigned int n_dofs() const override;

protected:
    const SmartPointer<const FiniteElement<dim>>  fe;
    const SmartPointer<const Quadrature<dim>>     quadrature;
    const SmartPointer<const Quadrature<dim - 1>> face_quadrature;
    DoFHandler<dim>                               dof_handler;
    Vector<double>                                solution;

    virtual void assemble_rhs(Vector<double> &rhs) const = 0;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs;

    virtual void apply_boundary_conditions() = 0;
    void setup_system();
    void assemble_system();
    void solve_system();
};


// CONSTRUCTOR
template <int dim>
Solver<dim>::Solver(Triangulation<dim> &       triangulation_,
                    const FiniteElement<dim> & fe_,
                    const Quadrature<dim> &    quadrature_,
                    const Quadrature<dim - 1> &face_quadrature_)
        : Base<dim>(triangulation_)
        , fe(&fe_)
        , quadrature(&quadrature_)
        , face_quadrature(&face_quadrature_)
        , dof_handler(triangulation_)
{}



#endif //ION_PROPULSION_SOLVER_H
