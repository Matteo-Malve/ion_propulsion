#include "Manifold.h"


template <int dim>
class Problem
{
public:
    Problem(Constants _constants);
    void run();

private:
    Constants constants;
    void create_mesh();

    void setup_primal_system();
    void setup_dual_system();
    void assemble_primal_system();
    void assemble_dual_system();
    void solve_primal();
    void solve_dual();
    void output_primal_results();
    void output_dual_results();

    double estimate_error(Vector<float> &error_indicators) const;
    void refine_mesh();

    Triangulation<dim> triangulation;

    FE_Q<dim> primal_fe;
    FE_Q<dim> dual_fe;

    QGauss<dim> dual_quadrature;

    DoFHandler<dim> dual_dof_handler;

    DoFHandler<dim> primal_dof_handler;

    AffineConstraints<double> primal_constraints; // to deal with hanging nodes
    SparsityPattern primal_sparsity_pattern;
    SparseMatrix<double> primal_system_matrix;

    Vector<double> primal_rhs;
    Vector<double> uh0;
    Vector<double> primal_solution;

    AffineConstraints<double> dual_constraints; // to deal with hanging nodes

    SparsityPattern dual_sparsity_pattern;
    SparseMatrix<double> dual_system_matrix;

    Vector<double> dual_solution;
    Vector<double> dual_rhs;

    RightHandSide<dim> rhs_function; // At the moment, this is equivalent to Functions::ZeroFunction<dim>()

    unsigned int cycle = 0;

    const unsigned int pre_refinement_steps = 0;
    const unsigned int max_refinements = 9;

    Timer timer;

    class Gradient;
    class IonizationArea;
};

// CONSTRUCTOR
template <int dim>
Problem<dim>::Problem(Constants _constants)
  : primal_fe(1)
  , dual_fe(2)
  , dual_quadrature(dual_fe.degree + 1)
  , primal_dof_handler(triangulation)
  , dual_dof_handler(triangulation)
  , constants(_constants)
  //, mapping()
{}

// INNER inhereted CLASSES


template <int dim>
class Problem<dim>::Gradient : public DataPostprocessorVector<dim>{
private:
    Constants constants;
public:
    Gradient(Constants _constants) : constants(_constants), DataPostprocessorVector<dim>("Electric_Field", update_gradients) {}

    virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                                       std::vector<Vector<double>> &computed_quantities) const override{

        AssertDimension(input_data.solution_gradients.size(), computed_quantities.size()); // size check

        for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
        {
            AssertDimension(computed_quantities[p].size(), dim); // dimension check
            for (unsigned int d = 0; d < dim; ++d)
                computed_quantities[p][d] = -input_data.solution_gradients[p][d];
        }
    }
};

template <int dim>
class Problem<dim>::IonizationArea : public DataPostprocessorScalar<dim>
{
private:
    Constants constants;
public:
    IonizationArea(Constants _constants): constants(_constants), DataPostprocessorScalar<dim>("normalized_overfield", update_gradients) {}

    virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                                       std::vector<Vector<double>> &computed_quantities) const override
    {
        auto E_ON = constants.E_ON;
        AssertDimension(input_data.solution_values.size(), computed_quantities.size()); // Size check

        for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
        {
            AssertDimension(computed_quantities[p].size(), 1);
            double magnitude = 0;
            for (unsigned int d = 0; d < dim; ++d)
                magnitude += input_data.solution_gradients[p][d] * input_data.solution_gradients[p][d];

            const double overfield = std::max(0., std::sqrt(magnitude) - E_ON);

            computed_quantities[p](0) = overfield / E_ON;
        }
    }
};



template <int dim>
class Evaluate_Rg : public Function<dim>
{
private:
    Constants constants;
public:
    Evaluate_Rg(Constants _constants): Function<dim>(), constants(_constants) {};
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override
    {
        auto X = constants.X;
        auto L = constants.L;
        auto R = constants.R;
        auto Vmax = constants.Vmax;
        (void)component;

        const auto y = p[1];
        const auto x = p[0];

        if (x <= (X - L / 2. - 2. * R) || x >= (X + L / 2. + 2. * R) || y >= 3. * R)
            return 0.;

        double Rg = 0.;
        double d = 0.;
        const double left_center = X - L / 2. + R;
        const double right_center = X + L / 2. - R;

        if (x >= right_center)
            d = sqrt((x - right_center) * (x - right_center) + y * y);
        else if (x < right_center && x > left_center)
            d = y;
        else if (x <= left_center)
            d = sqrt((x - left_center) * (x - left_center) + y * y);
        else
            cout << "ERROR! Not implemented!" << endl;

        d = std::max(d, R); /* As the boundary only approximates a circle,
                             * some points on the emitter edge will be inside the radius,
                             * and would produce an Rg higher than Vmax with the current function choice */

        if (d < 3. * R)
            Rg = Vmax * (1. - (d - R) / (2. * R)) * (1. - (d - R) / (2. * R)); //* (1. - (d-R)/(2.*R))

        if (Rg > Vmax)
            cout << "Error! Rg is " << Rg << " in " << x << ", " << y << endl;

        return Rg;
    }
};
