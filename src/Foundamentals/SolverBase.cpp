#include "SolverBase.h"

template <int dim>
void Solverbase<dim>::run()
{
    const unsigned int max_refinement = datafile("Mesh/Mesh_Refinement/max_refinement",20);
    const unsigned int min_refinement = datafile("Mesh/Mesh_Refinement/min_refinement",0);
    cout<<"\nSettings: "<<endl;
    cout<<" Min Refinemente = "<<min_refinement<<endl;
    cout << " Simulating: # refinements: "<<Nmax<< endl;

    while (cycle <= Nmax)
    {
        if (cycle == 0)
            create_mesh();
        else
            refine_grid(min_refinement, std::min(cycle + 1,max_refinement) );
        setup_system();

        // Forzante: imposta pari a 0
        assemble_rhs();

        // Permittività elettrica del vuoto:
        const double eps0 = 8.854*1e-12; // [F/m]
        // Matrice: matrice di Laplace, moltiplicata per la permittività dell'aria
        system_matrix.copy_from(laplace_matrix);
        system_matrix *= 1.0006*eps0*1e-3;

        constraints.condense(system_matrix, system_rhs);

        apply_boundary_conditions();


        std::cout << std::endl
                  << "Cycle " << cycle << ':' << std::endl
                  << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

        solve();
        output_results();

        std::cout << "   Elapsed CPU time: " << timer.cpu_time() << " seconds.\n";

        cycle++;
    }
}


template<int dim>
void Solverbase<dim>::apply_boundary_conditions() {
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    VectorTools::interpolate_boundary_values(dof_handler,
                                             1, // Boundary corrispondente all'emettitore, definito sopra
                                             Functions::ConstantFunction<dim>(2.e+4), // Valore di potenziale all'emettitore (20 kV)
                                             emitter_boundary_values);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             2,  // Boundary corrispondente al collettore, definito sopra
                                             Functions::ConstantFunction<dim>(1.6e+4), // Valore di potenziale al collettore (0 V)
                                             //DirichletBoundaryValuesDX<dim>(),
                                             collector_boundary_values);

    /* Le condizioni sopra sono condizioni di Dirichlet
    Su il restante bordo del dominio, dove non è imposta esplicitamente nessuna condizione,
    viene imposta automaticamente da deal.II una condizione di Neumann: gradiente nullo */

    MatrixTools::apply_boundary_values(emitter_boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);

    MatrixTools::apply_boundary_values(collector_boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}

template <int dim>
void Solverbase<dim>::create_mesh()
{
    CreateGrid(triangulation);

    // ----- Useless now that we fix boundary ids in gmsh
    //SetManifoldsAndBoundaries(triangulation, std::min(mesh_height,collector_height), electrode_distance, wire_radius);



    //ckeck_boundary_ids(triangulation);

    print_mesh_info(triangulation, "grid.svg");

    //triangulation.refine_global(1);

    //print_mesh_info(triangulation, "refined_grid.svg");

}

template <int dim>
void Solverbase<dim>::setup_system()
{
    dof_handler.distribute_dofs(fe);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
            /*keep_constrained_dofs = */ true);

    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);

    MatrixCreator::create_laplace_matrix(dof_handler,
                                         QGauss<dim>(fe.degree + 1),
                                         laplace_matrix);
}

template <int dim>
void Solverbase<dim>::solve()
{
    // Tolleranza e massimo numero di iterazioni
    unsigned int it_max = 1000;
    double tol = 1e-14;

    /*
     La tolleranza può anche essere definita in base alla norma l2
     del right hand side (che però nel nostro caso è assunto nullo):

     double tol = 1e-6 * system_rhs.l2_norm();

     */

    // Tipo di solver: CG (Conjugate Gradient)
    SolverControl            solver_control(it_max, tol);
    SolverCG<Vector<double>> solver(solver_control);

    // Parametro di rilassamento e precondizionatore
    // Sono anche disponibili Jacobi e parametrizzazione LU (vedi sotto)
    double relaxation_parameter = 1.2; // Questo parametro può variare nell'intervallo [1,2)
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, relaxation_parameter);

    /* Alternative:
     PreconditionJacobi<SparseMatrix<double>> preconditioner;
     preconditioner.initialize(system_matrix);

     oppure:

     SparseILU<SparseMatrix<double>> preconditioner;
     preconditioner.initialize(system_matrix);
     */

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    // Stampa a schermo del numero di iterazioni richieste
    std::cout << "   " << solver_control.last_step()
              << " CG iterations needed to obtain convergence." << std::endl;

    /* Il numero di iterazioni e il tempo CPU (sempre stampato a schermo)
     sono utili per definire la bontà di un solutore, in particolare
     studiando il loro andamento al variare degli N gradi di libertà
     della griglia per vedere se seguono un andamento proporzionale
     a N o N^2, per esempio
     */

    constraints.distribute(solution);
}

template <int dim>
void Solverbase<dim>::refine_grid(const unsigned int min_grid_level,
                               const unsigned int max_grid_level)
{
// Vettore errore (funziona solo con float, per quello che ho trovato)
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    // Stimatore dell'errore, studiato per l'equazione di Laplace
    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1), // Formula di quadratura di Gauss di grado pari agli elementi finiti + 1
                                       {},
                                       solution,
                                       estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.05 + 0.02*cycle,  // % celle con errore maggiore (da rifinire)
                                                    0.1); // % celle con errore minore (da "ri-unire")

    if (triangulation.n_levels() > max_grid_level)
        for (const auto &cell : triangulation.active_cell_iterators_on_level(max_grid_level))
            cell->clear_refine_flag();
    for (const auto &cell : triangulation.active_cell_iterators_on_level(min_grid_level))
        cell->clear_coarsen_flag();

    SolutionTransfer<dim> solution_trans(dof_handler);

    Vector<double> previous_solution;
    previous_solution = solution;
    triangulation.prepare_coarsening_and_refinement();
    solution_trans.prepare_for_coarsening_and_refinement(previous_solution);

    triangulation.execute_coarsening_and_refinement();
    setup_system();

    solution_trans.interpolate(previous_solution, solution);
    constraints.distribute(solution);
}

template <int dim>
void Solverbase<dim>::output_results()
{
    // Inizializzazione:
    if (cycle == 0) {
        values.reinit(Nmax+2);
        values(0) = 0;
    }

    /* Valore del potenziale in (1,0.2)
     (precedentemente usato come criterio di convergenza)

    double y = VectorTools::point_value(dof_handler, solution, {1.,0.2});
    std::cout << "   Potential at (1,0.2): " << y << std::endl;

    */

    /* Se, per qualche motivo, potesse essere utile valutare la
     convergenza sulla media del potenziale nel dominio,
     questa è la funzione per approssimarla:

     double x = VectorTools::compute_mean_value (dof_handler,QGauss<2>(fe.degree + 1),solution,0);

     */

    // !!! Criterio di convergenza
    Point<dim> evaluation_point(1.,0.2);

    // Print sample V
    Evaluation::PointValueEvaluation<dim> postprocessor(evaluation_point);
    double x_ = postprocessor(dof_handler,solution);
    std::cout << "   Potential at (" << evaluation_point[0] << "," << evaluation_point[1] << "): "
              << std::scientific << x_ << std::defaultfloat << std::endl;

    // Print sample E
    Tensor<1,dim>  E_ = VectorTools::point_gradient(dof_handler, solution, evaluation_point);
    double x = L2Norm(E_);
    std::cout << "   Field magnitude at (" << evaluation_point[0] << "," << evaluation_point[1] << "): " << x << std::endl;

    values(cycle+1) = x;

    // Se la condizione è verificata, si considera si sia raggiunta convergenza
    if ( std::fabs(values(cycle) - x) <= conv_tol*std::fabs(x))
        cycle = Nmax;

    // Scrittura su file della soluzione ad ogni ciclo:
    GradientPostprocessor<dim> gradient_postprocessor;
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "Potential");
    data_out.add_data_vector (solution, gradient_postprocessor);
    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
    data_out.write_vtu(output);


    // A convergenza raggiunta, stampa a schermo i risultati
    if (cycle == Nmax) {

        Point<dim> sample(wire_radius, 0.);
        Tensor<1,dim>  E = VectorTools::point_gradient(dof_handler, solution, sample);
        std::cout << "   Electric field in ("<< sample[0] << "," << sample[1] << "): " << -E << ", magnitude: " << L2Norm(E) << std::endl;
    }
}





// #######################################
// Template initialization
// #######################################
template class Solverbase<2>;

