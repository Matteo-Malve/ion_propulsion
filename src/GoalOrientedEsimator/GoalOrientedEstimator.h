#ifndef ION_PROPULSION_GOALORIENTEDESTIMATOR_H
#define ION_PROPULSION_GOALORIENTEDESTIMATOR_H
#include "../includes&parameters_setup.h"
#include "../Foundamentals/Problem.h"

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2002 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, ETH Zurich, 2002
 */


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <numeric>


#include "GOE_Framework.h" // Tutti gli altri a catena

/*

int main()
{
    try
    {
        using namespace Step14;

        const unsigned int                 dim = 2;
        Framework<dim>::ProblemDescription descriptor;

        descriptor.refinement_criterion =
                Framework<dim>::ProblemDescription::dual_weighted_error_estimator;

        descriptor.primal_fe_degree = 1;
        descriptor.dual_fe_degree   = 2;

        descriptor.data =
                std::make_unique<Data::SetUp<Data::Exercise_2_3<dim>, dim>>();

        const Point<dim> evaluation_point(0.75, 0.75);
        descriptor.dual_functional =
                std::make_unique<DualFunctional::PointValueEvaluation<dim>>(
                        evaluation_point);

        Evaluation::PointValueEvaluation<dim> postprocessor1(evaluation_point);
        Evaluation::GridOutput<dim>           postprocessor2("grid");

        descriptor.evaluator_list.push_back(&postprocessor1);
        descriptor.evaluator_list.push_back(&postprocessor2);

        descriptor.max_degrees_of_freedom = 20000;

        Framework<dim>::run(descriptor);
    }

    catch (std::exception &exc)
    {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }

    return 0;
}
 */

#endif //ION_PROPULSION_GOALORIENTEDESTIMATOR_H
