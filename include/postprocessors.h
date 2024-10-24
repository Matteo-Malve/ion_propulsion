#ifndef POSTPROCESSORS_H
#define POSTPROCESSORS_H

#include <deal.II/numerics/data_postprocessor.h>

using namespace dealii;


template <int dim>
class ElectricFieldPostprocessor : public DataPostprocessorVector<dim>
{
public:
  ElectricFieldPostprocessor ():  DataPostprocessorVector<dim> ("Electric_Field", update_gradients) {}

  virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
									 std::vector<Vector<double> > &computed_quantities) const override 
    {
    AssertDimension (input_data.solution_gradients.size(), computed_quantities.size()); // size check
    for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p){
      AssertDimension (computed_quantities[p].size(), dim); // dimension check
      for (unsigned int d=0; d<dim; ++d)
        computed_quantities[p][d] = -input_data.solution_gradients[p][d];
      }
    }
};


template <int dim>
class HomogeneousFieldPostprocessor : public DataPostprocessorVector<dim>
{
public:
  HomogeneousFieldPostprocessor ():  DataPostprocessorVector<dim> ("minus_grad_uh0", update_gradients) {}

  virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
									 std::vector<Vector<double> > &computed_quantities) const override 
    {
    AssertDimension (input_data.solution_gradients.size(), computed_quantities.size()); // size check
    for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p){
      AssertDimension (computed_quantities[p].size(), dim); // dimension check
      for (unsigned int d=0; d<dim; ++d)
        computed_quantities[p][d] = -input_data.solution_gradients[p][d];
      }
    }
};




#endif