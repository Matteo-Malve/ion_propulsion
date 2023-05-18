#ifndef GETPOT_GOE_LAPLACESOLVER_H
#define GETPOT_GOE_LAPLACESOLVER_H

#include "GOE_Evaluation.h"

namespace GOE {
    using namespace dealii;
    namespace LaplaceSolver {

        template<int dim>
        class Base {
        public:
            Base(Triangulation <dim> &coarse_grid);

            virtual ~Base() = default;

            virtual void solve_problem() = 0;

            virtual void postprocess(
                    const Evaluation::EvaluationBase <dim> &postprocessor) const = 0;

            virtual void refine_grid() = 0;

            virtual unsigned int n_dofs() const = 0;

            virtual void set_refinement_cycle(const unsigned int cycle);

            virtual void output_solution() const = 0;

        protected:
            const SmartPointer <Triangulation<dim>> triangulation;

            unsigned int refinement_cycle;
        };


        template<int dim>
        Base<dim>::Base(Triangulation <dim> &coarse_grid)
                : triangulation(&coarse_grid), refinement_cycle(numbers::invalid_unsigned_int) {}


        template<int dim>
        void Base<dim>::set_refinement_cycle(const unsigned int cycle) {
            refinement_cycle = cycle;
        }


        template<int dim>
        class Solver : public virtual Base<dim> {
        public:
            Solver(Triangulation <dim> &triangulation,
                   const FiniteElement <dim> &fe,
                   const Quadrature <dim> &quadrature,
                   const Quadrature<dim - 1> &face_quadrature,
                   const Function <dim> &boundary_values);

            virtual ~Solver() override;

            virtual void solve_problem() override;

            virtual void postprocess(
                    const Evaluation::EvaluationBase <dim> &postprocessor) const override;

            virtual unsigned int n_dofs() const override;

        protected:
            const SmartPointer<const FiniteElement <dim>> fe;
            const SmartPointer<const Quadrature <dim>> quadrature;
            const SmartPointer< const Quadrature<dim - 1>>
            face_quadrature;
            DoFHandler <dim> dof_handler;
            Vector<double> solution;
            const SmartPointer<const Function <dim>> boundary_values;

            virtual void assemble_rhs(Vector<double> &rhs) const = 0;

        private:
            struct LinearSystem {
                LinearSystem(const DoFHandler <dim> &dof_handler);

                void solve(Vector<double> &solution) const;

                AffineConstraints<double> hanging_node_constraints;
                SparsityPattern sparsity_pattern;
                SparseMatrix<double> matrix;
                Vector<double> rhs;
            };


            struct AssemblyScratchData {
                AssemblyScratchData(const FiniteElement <dim> &fe,
                                    const Quadrature <dim> &quadrature);

                AssemblyScratchData(const AssemblyScratchData &scratch_data);

                FEValues <dim> fe_values;
            };

            struct AssemblyCopyData {
                FullMatrix<double> cell_matrix;
                std::vector <types::global_dof_index> local_dof_indices;
            };


            void assemble_linear_system(LinearSystem &linear_system);

            void local_assemble_matrix(
                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                    AssemblyScratchData &scratch_data,
                    AssemblyCopyData &copy_data) const;


            void copy_local_to_global(const AssemblyCopyData &copy_data,
                                      LinearSystem &linear_system) const;
        };


        template<int dim>
        Solver<dim>::Solver(Triangulation <dim> &triangulation,
                            const FiniteElement <dim> &fe,
                            const Quadrature <dim> &quadrature,
                            const Quadrature<dim - 1> &face_quadrature,
                            const Function <dim> &boundary_values)
                : Base<dim>(triangulation), fe(&fe), quadrature(&quadrature), face_quadrature(&face_quadrature),
                  dof_handler(triangulation), boundary_values(&boundary_values) {}


        template<int dim>
        Solver<dim>::~Solver() {
            dof_handler.clear();
        }


        template<int dim>
        void Solver<dim>::solve_problem() {
            dof_handler.distribute_dofs(*fe);
            solution.reinit(dof_handler.n_dofs());

            LinearSystem linear_system(dof_handler);
            assemble_linear_system(linear_system);
            linear_system.solve(solution);
        }


        template<int dim>
        void Solver<dim>::postprocess(
                const Evaluation::EvaluationBase <dim> &postprocessor) const {
            postprocessor(dof_handler, solution);
        }


        template<int dim>
        unsigned int Solver<dim>::n_dofs() const {
            return dof_handler.n_dofs();
        }


        template<int dim>
        void Solver<dim>::assemble_linear_system(LinearSystem &linear_system) {
            Threads::Task<void> rhs_task =
                    Threads::new_task(&Solver<dim>::assemble_rhs, *this, linear_system.rhs);

            auto worker =
                    [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
                           AssemblyScratchData &scratch_data,
                           AssemblyCopyData &copy_data) {
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


        template<int dim>
        Solver<dim>::AssemblyScratchData::AssemblyScratchData(
                const FiniteElement <dim> &fe,
                const Quadrature <dim> &quadrature)
                : fe_values(fe, quadrature, update_gradients | update_JxW_values) {}


        template<int dim>
        Solver<dim>::AssemblyScratchData::AssemblyScratchData(
                const AssemblyScratchData &scratch_data)
                : fe_values(scratch_data.fe_values.get_fe(),
                            scratch_data.fe_values.get_quadrature(),
                            update_gradients | update_JxW_values) {}


        template<int dim>
        void Solver<dim>::local_assemble_matrix(
                const typename DoFHandler<dim>::active_cell_iterator &cell,
                AssemblyScratchData &scratch_data,
                AssemblyCopyData &copy_data) const {
            const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
            const unsigned int n_q_points = quadrature->size();

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


        template<int dim>
        void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data,
                                               LinearSystem &linear_system) const {
            for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i)
                for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j)
                    linear_system.matrix.add(copy_data.local_dof_indices[i],
                                             copy_data.local_dof_indices[j],
                                             copy_data.cell_matrix(i, j));
        }


        template<int dim>
        Solver<dim>::LinearSystem::LinearSystem(const DoFHandler <dim> &dof_handler) {
            hanging_node_constraints.clear();

            void (*mhnc_p)(const DoFHandler <dim> &, AffineConstraints<double> &) =
            &DoFTools::make_hanging_node_constraints;

            Threads::Task<void> side_task =
                    Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints);

            DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
            DoFTools::make_sparsity_pattern(dof_handler, dsp);


            side_task.join();

            hanging_node_constraints.close();
            hanging_node_constraints.condense(dsp);
            sparsity_pattern.copy_from(dsp);

            matrix.reinit(sparsity_pattern);
            rhs.reinit(dof_handler.n_dofs());
        }


        template<int dim>
        void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const {
            SolverControl solver_control(5000, 1e-12);
            SolverCG <Vector<double>> cg(solver_control);

            PreconditionSSOR <SparseMatrix<double>> preconditioner;
            preconditioner.initialize(matrix, 1.2);

            cg.solve(matrix, solution, rhs, preconditioner);

            hanging_node_constraints.distribute(solution);
        }


        template<int dim>
        class PrimalSolver : public Solver<dim> {
        public:
            PrimalSolver(Triangulation <dim> &triangulation,
                         const FiniteElement <dim> &fe,
                         const Quadrature <dim> &quadrature,
                         const Quadrature<dim - 1> &face_quadrature,
                         const Function <dim> &rhs_function,
                         const Function <dim> &boundary_values);

            virtual void output_solution() const override;

        protected:
            const SmartPointer<const Function <dim>> rhs_function;

            virtual void assemble_rhs(Vector<double> &rhs) const override;
        };


        template<int dim>
        PrimalSolver<dim>::PrimalSolver(Triangulation <dim> &triangulation,
                                        const FiniteElement <dim> &fe,
                                        const Quadrature <dim> &quadrature,
                                        const Quadrature<dim - 1> &face_quadrature,
                                        const Function <dim> &rhs_function,
                                        const Function <dim> &boundary_values)
                : Base<dim>(triangulation), Solver<dim>(triangulation,
                                                        fe,
                                                        quadrature,
                                                        face_quadrature,
                                                        boundary_values), rhs_function(&rhs_function) {}


        template<int dim>
        void PrimalSolver<dim>::output_solution() const {
            DataOut <dim> data_out;
            data_out.attach_dof_handler(this->dof_handler);
            data_out.add_data_vector(this->solution, "solution");
            data_out.build_patches();

            std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
                              ".vtu");
            data_out.write(out, DataOutBase::vtu);
        }


        template<int dim>
        void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const {
            FEValues <dim> fe_values(*this->fe,
                                     *this->quadrature,
                                     update_values | update_quadrature_points |
                                     update_JxW_values);

            const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
            const unsigned int n_q_points = this->quadrature->size();

            Vector<double> cell_rhs(dofs_per_cell);
            std::vector<double> rhs_values(n_q_points);
            std::vector <types::global_dof_index> local_dof_indices(dofs_per_cell);

            for (const auto &cell: this->dof_handler.active_cell_iterators()) {
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


        template<int dim>
        class RefinementGlobal : public PrimalSolver<dim> {
        public:
            RefinementGlobal(Triangulation <dim> &coarse_grid,
                             const FiniteElement <dim> &fe,
                             const Quadrature <dim> &quadrature,
                             const Quadrature<dim - 1> &face_quadrature,
                             const Function <dim> &rhs_function,
                             const Function <dim> &boundary_values);

            virtual void refine_grid() override;
        };


        template<int dim>
        RefinementGlobal<dim>::RefinementGlobal(
                Triangulation <dim> &coarse_grid,
                const FiniteElement <dim> &fe,
                const Quadrature <dim> &quadrature,
                const Quadrature<dim - 1> &face_quadrature,
                const Function <dim> &rhs_function,
                const Function <dim> &boundary_values)
                : Base<dim>(coarse_grid), PrimalSolver<dim>(coarse_grid,
                                                            fe,
                                                            quadrature,
                                                            face_quadrature,
                                                            rhs_function,
                                                            boundary_values) {}


        template<int dim>
        void RefinementGlobal<dim>::refine_grid() {
            this->triangulation->refine_global(1);
        }


        template<int dim>
        class RefinementKelly : public PrimalSolver<dim> {
        public:
            RefinementKelly(Triangulation <dim> &coarse_grid,
                            const FiniteElement <dim> &fe,
                            const Quadrature <dim> &quadrature,
                            const Quadrature<dim - 1> &face_quadrature,
                            const Function <dim> &rhs_function,
                            const Function <dim> &boundary_values);

            virtual void refine_grid() override;
        };


        template<int dim>
        RefinementKelly<dim>::RefinementKelly(
                Triangulation <dim> &coarse_grid,
                const FiniteElement <dim> &fe,
                const Quadrature <dim> &quadrature,
                const Quadrature<dim - 1> &face_quadrature,
                const Function <dim> &rhs_function,
                const Function <dim> &boundary_values)
                : Base<dim>(coarse_grid), PrimalSolver<dim>(coarse_grid,
                                                            fe,
                                                            quadrature,
                                                            face_quadrature,
                                                            rhs_function,
                                                            boundary_values) {}


        template<int dim>
        void RefinementKelly<dim>::refine_grid() {
            Vector<float> estimated_error_per_cell(
                    this->triangulation->n_active_cells());
            KellyErrorEstimator<dim>::estimate(
                    this->dof_handler,
                    QGauss<dim - 1>(this->fe->degree + 1),
                    std::map<types::boundary_id, const Function <dim> *>(),
                    this->solution,
                    estimated_error_per_cell);
            GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
                                                            estimated_error_per_cell,
                                                            0.3,
                                                            0.03);
            this->triangulation->execute_coarsening_and_refinement();
        }


        template<int dim>
        class RefinementWeightedKelly : public PrimalSolver<dim> {
        public:
            RefinementWeightedKelly(Triangulation <dim> &coarse_grid,
                                    const FiniteElement <dim> &fe,
                                    const Quadrature <dim> &quadrature,
                                    const Quadrature<dim - 1> &face_quadrature,
                                    const Function <dim> &rhs_function,
                                    const Function <dim> &boundary_values,
                                    const Function <dim> &weighting_function);

            virtual void refine_grid() override;

        private:
            const SmartPointer<const Function <dim>> weighting_function;
        };


        template<int dim>
        RefinementWeightedKelly<dim>::RefinementWeightedKelly(
                Triangulation <dim> &coarse_grid,
                const FiniteElement <dim> &fe,
                const Quadrature <dim> &quadrature,
                const Quadrature<dim - 1> &face_quadrature,
                const Function <dim> &rhs_function,
                const Function <dim> &boundary_values,
                const Function <dim> &weighting_function)
                : Base<dim>(coarse_grid), PrimalSolver<dim>(coarse_grid,
                                                            fe,
                                                            quadrature,
                                                            face_quadrature,
                                                            rhs_function,
                                                            boundary_values), weighting_function(&weighting_function) {}


        template<int dim>
        void RefinementWeightedKelly<dim>::refine_grid() {
            Vector<float> estimated_error_per_cell(
                    this->triangulation->n_active_cells());
            std::map<types::boundary_id, const Function <dim> *> dummy_function_map;
            KellyErrorEstimator<dim>::estimate(this->dof_handler,
                                               *this->face_quadrature,
                                               dummy_function_map,
                                               this->solution,
                                               estimated_error_per_cell);

            for (const auto &cell: this->dof_handler.active_cell_iterators())
                estimated_error_per_cell(cell->active_cell_index()) *=
                        weighting_function->value(cell->center());

            GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
                                                            estimated_error_per_cell,
                                                            0.3,
                                                            0.03);
            this->triangulation->execute_coarsening_and_refinement();
        }

    } // namespace LaplaceSolver
}

#endif //GETPOT_GOE_LAPLACESOLVER_H
