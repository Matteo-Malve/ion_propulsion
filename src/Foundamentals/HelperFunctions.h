#ifndef MESH_CC_HELPERFUCTIONS_H
#define MESH_CC_HELPERFUCTIONS_H

#include "../includes&parameters_setup.h"

// Funzione per stampare a schermo alcune caratteristiche della griglia
template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename);

// Funzione per calcolo norma di tensore
template <int dim>
double L2Norm(const Tensor<1,dim> &input);

// Funzione usata in post-processing per estrarre il campo elettrico come -gradiente del potenziale
template <int dim>
class GradientPostprocessor : public DataPostprocessorVector<dim>
{
public:
    GradientPostprocessor ():  DataPostprocessorVector<dim> ("E", update_gradients) {}

    virtual void evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                                       std::vector<Vector<double> > &computed_quantities) const override {
        AssertDimension (input_data.solution_gradients.size(), computed_quantities.size()); // size check

        for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
        {
            AssertDimension (computed_quantities[p].size(), dim); // dimension check
            for (unsigned int d=0; d<dim; ++d)
                computed_quantities[p][d] = -input_data.solution_gradients[p][d];
        }
    }
};

template <int dim>
void ckeck_boundary_ids(const Triangulation<dim> &triangulation);

#endif //MESH_CC_HELPERFUCTIONS_H
