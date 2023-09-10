#ifndef ION_PROPULSION_CC_HELPERFUCTIONS_H
#define ION_PROPULSION_CC_HELPERFUCTIONS_H

#include "../includes&parameters_setup.h"
static GetPot redefined_4_datafile("../data_setup");

static auto evaluate_grad_Rg = [](double x, double y) {
    double r = sqrt(x * x + y * y);
    Tensor<1,2> grad_Rg;
    double Ve = 20000;
    double Re = 250e-6;

    // Gradient computed analytically by hand. Check Evaluate_Rg.h for the primitive function.
    grad_Rg[0] = 0;
    grad_Rg[1] = 0;
    if(r<2*Re) {
        grad_Rg[0] = -Ve / Re * x / r;
        grad_Rg[1] = -Ve / Re * y / r;
    }
    return grad_Rg;
};


// To compute L2-norm of a tensor
template <int dim>
double L2Norm(const Tensor<1,dim> &input);

// To be used in post-processing to compute the electrical field as -gradient of the potential
template <int dim>
class GradientPostprocessor : public DataPostprocessorVector<dim>
{
public:
    GradientPostprocessor ():  DataPostprocessorVector<dim> ("E", update_gradients) {}

    virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                                       std::vector<Vector<double> > &computed_quantities) const override {
        AssertDimension (input_data.solution_gradients.size(), computed_quantities.size()); // Size check

        for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
        {
            AssertDimension (computed_quantities[p].size(), dim); // Dimensionality check
            for (unsigned int d=0; d<dim; ++d)
                computed_quantities[p][d] = -input_data.solution_gradients[p][d];
        }
    }
};

// Our target function, retrieves the potential area of ionized air, basing on a threshold
template <int dim>
void ionization_area(const Triangulation<dim> &triangulation, const DoFHandler<dim> &dof_handler, const Vector<double> &solution);


/* Helper function, gives errors on less recent versions of deal.ii
 * It outputs a .vtu with a field dedicated to the boundary ids.
 * It is helpful to visualize if you've correctly set the desired BCs.

template <int dim>
void check_boundary_ids(const Triangulation<dim> &triangulation);

*/

#endif //ION_PROPULSION_CC_HELPERFUCTIONS_H
