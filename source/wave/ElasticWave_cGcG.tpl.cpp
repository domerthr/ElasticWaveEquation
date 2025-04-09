// PROJECT includes
#include <wave/ElasticWave_cGcG.tpl.hh>

#include <wave/InitialValue/InitialValue_Selector.tpl.hh>
#include <wave/FinalValue/FinalValue_Selector.tpl.hh>
#include <wave/Force/Force_Selector.tpl.hh>
#include <wave/Boundary/BoundaryValue_Selector.tpl.hh>
#include <wave/ExactSolution/ExactSolution_Selector.tpl.hh>

#include <fest/fe/fe.tpl.hh>



// DEAL.II includes



// C++ includes



template <int dim>
ElasticWave_cGcG<dim>::ElasticWave_cGcG(unsigned int s, unsigned int r)
    :   grid(),
        dof_handler(grid),
        fe(std::make_shared< dealii::FESystem<dim> >(/*u*/ dealii::FE_Q<dim>(s), dim,
                                                     /*v*/ dealii::FE_Q<dim>(s), dim), 
           r)
{}

template <int dim>
inline dealii::SymmetricTensor<2, dim> ElasticWave_cGcG<dim>::get_strain(const dealii::Tensor<2,dim> &grad) {
    Assert(grad.dimension == dim, dealii::ExcInternalError());
  
    dealii::SymmetricTensor<2, dim> strain;
    for (unsigned int i = 0; i < dim; ++i)
      strain[i][i] = grad[i][i];

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i + 1; j < dim; ++j)
        strain[i][j] = (grad[i][j] + grad[j][i]) / 2;

    return strain;
}


template <int dim>
void ElasticWave_cGcG<dim>::set_input_parameters(
    std::shared_ptr< dealii::ParameterHandler > parameter_handler) {
    
    parameter_set = std::make_shared< wave::ParameterSet > (
		parameter_handler); 
}


template <int dim>
void ElasticWave_cGcG<dim>::init_functions() {

    // inital value function
    {
        wave::initial_value::Selector<dim> selector;
        selector.create_function(
                parameter_set->example,
                function.inital_value_function
        );
    }
    
    // final value function
    {
        wave::final_value::Selector<dim> selector;
        selector.create_function(
                parameter_set->example,
                function.final_value_function,
                parameter_set->T
        );
    }
    // force function
    {
        wave::force::Selector<dim> selector;
        selector.create_function(
                parameter_set->example,
                function.rhs_function
        );
    }
    // boundary value function
    {
        wave::boundary_value::Selector<dim> selector;
        selector.create_function(
                parameter_set->example,
                function.boundary_function
        );
    }
    // exact solution function
    {
        wave::exact_solution::Selector<dim> selector;
        selector.create_function(
                parameter_set->example,
                function.exact_solution
        );
    }
}

template <int dim>
void ElasticWave_cGcG<dim>::init_grid() {


    // generator to generate spatial triangulation
    wave::grid::TriaGenerator<dim> tria_gen;
    tria_gen.generate(grid.spatial(), parameter_set->example);

    // generate a Space-Time Triangulation
    grid.generate(parameter_set->t0, parameter_set->T);
    // refine the grid in (space, time)
    grid.refine_global(parameter_set->global_refinement_space, parameter_set->global_refinement_time);
    
}


template <int dim>
void ElasticWave_cGcG<dim>::setup_dofs() {


    // Distribute spatial and temporal dofs
    dof_handler.distribute_dofs(fe);


    std::vector<unsigned int> block_component(dim + dim, 0);
    std::fill(block_component.begin() + dim, block_component.end(), 1);


    // Renumbering spartial DoFs into displacement and velocity
    dealii::DoFRenumbering::component_wise(*(dof_handler.spatial()), block_component);
       
    // two blocks: displacement and velocity
    const std::vector<dealii::types::global_dof_index> dofs_per_block = 
          dealii::DoFTools::count_dofs_per_fe_block(*(dof_handler.spatial()), block_component);

    n_space_u = dofs_per_block[0];
    n_space_v = dofs_per_block[1];

    const unsigned int n_dofs_space = dof_handler.n_dofs_space();
    const unsigned int n_dofs_time  = dof_handler.n_dofs_time();

    // computing Initial Values
    initial_value.reinit(n_dofs_space);
    dealii::VectorTools::interpolate(*(dof_handler.spatial()), *function.inital_value_function, initial_value, dealii::ComponentMask());
    
    // set homogenouse Dirichlet boundary constraints
    constraints = std::make_shared< dealii::AffineConstraints<double> > ();
    fest::cgcg::VectorTools::interpolate_boundary_values(dof_handler, 0,*function.boundary_function, constraints);

    constraints->close();


    // Construct the space-time sparsity pattern for this slab
    dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs_spacetime());
    fest::cgcg::DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    // print sparsity_pattern of the space-time slab
    // // std::ofstream out_sparsity("../sparsity_pattern.svg");
    // // sparsity_pattern.print_svg(out_sparsity);

    // reinit the linear system
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs_spacetime());
    system_rhs.reinit(dof_handler.n_dofs_spacetime());

    
    refinement_table.add_value("dofs (space)", n_dofs_space);
    refinement_table.add_value("dofs (time)", n_dofs_time);


}


template <int dim>
void ElasticWave_cGcG<dim>::assemble_system() {
    
    // reset slab system matrix and right hand side
    system_matrix = 0;
    system_rhs = 0;

    // space - time quadratur and FEValues
    fest::QGauss<dim> quad(fe.spatial()->degree + 2, fe.temporal()->degree + 2);
    fest::cgcg::FEValues<dim> fe_values(fe, quad,
                                        dealii::update_values | dealii::update_gradients | 
                                        dealii::update_quadrature_points | dealii::update_JxW_values);


    const unsigned int dofs_per_spacetime_cell = fe.dofs_per_cell;

    // Number of temporal elements and index of the current elelemt for offset calculations
    unsigned int N = grid.temporal()->n_global_active_cells();
    unsigned int n;

    // space-time local matrix and rhs containing all temporal information but only one spatial element
    dealii::FullMatrix<double> local_matrix(N * dofs_per_spacetime_cell, 
                                            N * dofs_per_spacetime_cell);
    dealii::Vector<double> local_rhs(N * dofs_per_spacetime_cell);
    std::vector<dealii::types::global_dof_index> local_dof_indices(N * 
                                                                   dofs_per_spacetime_cell);

    // number of space-time and space quadratur points
    unsigned int n_quad_spacetime = fe_values.n_quadature_points;
    
    double lame_coefficient_lambda = (2 * parameter_set->poisson_ratio_nu * parameter_set->lame_coefficient_mu) / (1.0 - 2 * parameter_set->poisson_ratio_nu);
    
    dealii::SymmetricTensor<2, dim> identity = dealii::unit_symmetric_tensor<dim, double>();
    
    // First iterate over active spatial elements
    for (const auto &space_cell : dof_handler.spatial()->active_cell_iterators()) {

        // reset local contributions
        local_matrix = 0;
        local_rhs = 0;

        // recalculate local information for the current spatial element
        fe_values.reinit_space(space_cell);
        
        
        // iterate over all active temporal elements
        for (const auto &time_cell : dof_handler.temporal()->active_cell_iterators()) {
            n = time_cell->index();
            
            // recalculate the local information of the current temporla element
            fe_values.reinit_time(time_cell);

            // Get local space-time dof-indices for the current tensor product cell
            fe_values.get_local_dof_indices(local_dof_indices);

            
            // Iterate over all quadratur points except for jump terms
            for (unsigned int q = 0; q < n_quad_spacetime; ++q) {

                // set the time of the right hand side of the current quadratur point and get its location
                const double t_qq = fe_values.time_quadratur_point(q);
                function.rhs_function->set_time(t_qq);
                const auto &x_q = fe_values.space_quadratur_point(q);

                // Iterate over all space-time dofs of the current element
                for (unsigned int i = 0; i < dofs_per_spacetime_cell; ++i) {

                    local_rhs(i + n * dofs_per_spacetime_cell) += 
                        (function.rhs_function->value(x_q) *                               // f(t_qq, x_q)
                         fe_values.vector_value(displacement, i, q) *                      // φ^u_{i_time, i_space}(t_qq, x_q)
                         fe_values.JxW(q));                                                // dx dt
                    

                    for (unsigned int j = 0; j < dofs_per_spacetime_cell; ++j) {

                        local_matrix(j + n * dofs_per_spacetime_cell,
                                     i + n * dofs_per_spacetime_cell) += 
                            
                             (parameter_set->rho *                                          // ρ
                              fe_values.vector_dt(velocity, i, q) *                         // ∂_t φ^v_{i_time, i_space}(t_qq, x_q)
                              fe_values.vector_value(displacement, j, q)                    // φ^u_{j_time, j_space}(t_qq, x_q)
                              

                            + dealii::scalar_product(2 * parameter_set->lame_coefficient_mu *   // 2 * μ
                              get_strain(fe_values.vector_space_grad(displacement, i, q)),      // 1/2 * (∇_x φ^u_{j_time, j_space}(t_qq, x_q) + ∇_x^T φ^u_{j_time, j_space}(t_qq, x_q))
                              fe_values.vector_space_grad(displacement, j, q))                  // ∇_x φ^u_{i_time, i_space}(t_qq, x_q)
                              
                            
                            + dealii::scalar_product(lame_coefficient_lambda *                                          // λ
                              dealii::trace(get_strain(fe_values.vector_space_grad(displacement, i, q))) * identity,    // tr(1/2 * (∇_x φ^u_{j_time, j_space}(t_qq, x_q) + ∇_x^T φ^u_{j_time, j_space}(t_qq, x_q)))
                              fe_values.vector_space_grad(displacement, j, q))                                          // ∇_x φ^u_{i_time, i_space}(t_qq, x_q)
                              
                                                          
                            + fe_values.vector_dt(displacement, i, q) *                     // ∂_t φ^u_{i_time, i_space}(t_qq, x_q)
                              fe_values.vector_value(velocity, j, q)                        // φ^v_{j_time, j_space}(t_qq, x_q)
                                 

                            - fe_values.vector_value(velocity, i, q) *                      // φ^v_{i_time, i_space}(t_qq, x_q)
                              fe_values.vector_value(velocity, j, q))                       // φ^v_{j_time, j_space}(t_qq, x_q)

                            * fe_values.JxW(q);                                             // dx dt
                    } // dof j
                } // dof i
                    
            } // q


        } // cell time

        
        constraints->distribute_local_to_global(local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);


    } // cell space
    
    apply_temporal_conditions();
    

}

template <int dim>
void ElasticWave_cGcG<dim>::apply_temporal_conditions()
{
    dealii::AffineConstraints<double> constraints;
    constraints.close();


    if (parameter_set->junker) {    

        // apply initial value v(t_0) on linear system
        // BUT do nothing for u(t_0)
        apply_initial_values(n_space_u, dof_handler.n_dofs_space());
    
        // apply final value v(T) on linear system
        // BUT do nothing for u(T)
        apply_final_values();

    }
    else {

        // apply commen initial values u(t_0) and v(t_0) on linear system
        apply_initial_values(0, dof_handler.n_dofs_space());
        
    }   

}

template <int dim>
void ElasticWave_cGcG<dim>::apply_initial_values(unsigned int start, unsigned int end) {

    // clear from start to end all rows of system matrix
    for (unsigned int i = start; i < end; ++i){
        // iterate over all entries and set A_ij = δ_ij
        for (typename dealii::SparseMatrix<double>::iterator j = system_matrix.begin(i); j != system_matrix.end(i); ++j) {
            j->value() = 0.;
            system_matrix.set(i, i, 1.);
            system_rhs(i) = initial_value(i);
        }
   }

}

template <int dim>
void ElasticWave_cGcG<dim>::apply_final_values() {


    // computing Final Values
    dealii::Vector<double> final_value(dof_handler.n_dofs_space());
    dealii::VectorTools::interpolate(*(dof_handler.spatial()), *function.final_value_function, final_value, dealii::ComponentMask());

    // claer last N_space rows of system matrix and skip u-DoFs
    for (unsigned int i = dof_handler.n_dofs_spacetime() - n_space_v; i < dof_handler.n_dofs_spacetime(); ++i) {
        // iterate over all entries and set A_ij = δ_ij
        for (typename dealii::SparseMatrix<double>::iterator j = system_matrix.begin(i); j != system_matrix.end(i); ++j) {
            j->value() = 0.;
            system_matrix.set(i, i, 1.);
            system_rhs(i) = final_value(i + n_space_u + n_space_v - dof_handler.n_dofs_spacetime());
        }
    }

}



template <int dim>
void ElasticWave_cGcG<dim>::solve() {


    dealii::SparseDirectUMFPACK solver;
    solver.factorize(system_matrix);
    solver.vmult(solution, system_rhs);
    // dealii::SolverControl solver_control(10000, 1e-6, log_result=false);
    // //SolverFGMRES<BlockVector<double>> solver(solver_control);
    // //SolverGMRES<BlockVector<double>> solver(solver_control);
    // dealii::SolverBicgstab<dealii::Vector<double>> solver(solver_control);

    // dealii::PreconditionJacobi<dealii::SparseMatrix<double>> preconditioner;
    // preconditioner.initialize(system_matrix);

    // solver.solve(system_matrix, solution, system_rhs, preconditioner);


    constraints->distribute(solution);

}

template <int dim>
void ElasticWave_cGcG<dim>::output_results(const unsigned int refinement_cycle) {
    
    const std::string output_dir = "../output/cG-cG/example-" + std::to_string(parameter_set->example) + "/"
                                    + "final_condition=" + std::to_string(parameter_set->junker) + "/" 
                                    + "cylce=" + std::to_string(refinement_cycle) + "/";

    std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> 
            data_component_interpretation(2*dim, dealii::DataComponentInterpretation::component_is_part_of_vector);
    std::vector<std::string> solution_names(dim, "displacement");
    std::vector<std::string> solution_names_v(dim, "velocity");
    solution_names.insert(solution_names.end(), solution_names_v.begin(), solution_names_v.end());
    
    auto n_dofs = dof_handler.n_dofs_time();

    // time
    std::vector<dealii::Point<1>> time_support_points(dof_handler.temporal()->n_dofs());
    dealii::DoFTools::map_dofs_to_support_points(
        dealii::MappingQ1<1, 1>(), *dof_handler.temporal(), time_support_points);

    
    for (unsigned int i = 0; i < n_dofs; i++) {

        dealii::Point<1> point = time_support_points[i];
        
        dealii::DataOut<dim> data_out;
        data_out.attach_dof_handler(*dof_handler.spatial());
        dealii::Vector<double> local_solution;
        local_solution.reinit(dof_handler.n_dofs_space());

        fest::cgcg::VectorTools::extract_subvector_at_time_dof(solution, local_solution, i);
     
        data_out.add_data_vector(local_solution, solution_names, dealii::DataOut<dim>::type_dof_data, data_component_interpretation);
        data_out.build_patches();
        data_out.set_flags(dealii::DataOutBase::VtkFlags(point[0]));
        
        std::ostringstream filename;
        filename << "solution(" << fe.temporal()->degree << ")_t_" << i << ".vtk";

        std::ofstream output(output_dir + filename.str());
        data_out.write_vtk(output);

        output.close();

    }    
    
}

template <int dim>
void ElasticWave_cGcG<dim>::process_solution(unsigned int cycle) {


    fest::QGauss<dim> quad(fe.spatial()->degree + 2, fe.temporal()->degree + 2);
    L2_error = fest::cgcg::VectorTools::calculate_L2L2_squared_error<dim>(dof_handler, solution, *function.exact_solution, quad);


    
    L2_error = std::sqrt(L2_error);
    L2_error_vals.push_back(L2_error);

    const unsigned int n_dofs_space = dof_handler.n_dofs_space();
    const unsigned int n_dofs_time  = dof_handler.n_dofs_time();
    const unsigned int n_dofs = n_dofs_space * n_dofs_time;

    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", grid.n_active_cells());
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("dofs(space)", n_dofs_space);
    convergence_table.add_value("dofs(time)", n_dofs_time);
    convergence_table.add_value("L2", L2_error);
    
    
}





template <int dim>
void ElasticWave_cGcG<dim>::plot_solution(const unsigned int refinement_cycle) {

    const std::string output_dir = "../output/cG-cG/example-" + std::to_string(parameter_set->example) + "/"
                                    + "final_condition=" + std::to_string(parameter_set->junker) + "/" 
                                    + "cylce=" + std::to_string(refinement_cycle) + "/";

    std::string py_cmd2 = "python3 ../python/plot_solution.py " + output_dir + " " + std::to_string(parameter_set->example);
    system(py_cmd2.c_str());
}

template <int dim>
void ElasticWave_cGcG<dim>::print_convergence_table() {

    const std::string output_dir = "../output/cG-cG/example-" + std::to_string(parameter_set->example) + "/"
                                    + "final_condition=" + std::to_string(parameter_set->junker) + "/";

    convergence_table.set_precision("L2", 16);
    convergence_table.set_scientific("L2", true);
    convergence_table.set_tex_caption("cells", "\\# cells");
    convergence_table.set_tex_caption("dofs", "\\# dofs");
    convergence_table.set_tex_caption("dofs(space)", "\\# dofs space");
    convergence_table.set_tex_caption("dofs(time)", "\\# dofs time");
    convergence_table.set_tex_caption("L2", "@f$L^2@f$-error");
    convergence_table.set_tex_format("cells", "r");
    convergence_table.set_tex_format("dofs", "r");
    convergence_table.set_tex_format("dofs(space)", "r");
    convergence_table.set_tex_format("dofs(time)", "r");
    std::cout << std::endl;
    convergence_table.write_text(std::cout, dealii::TableHandler::TextOutputFormat::org_mode_table);

    std::ofstream output_1(output_dir + "convergence_table.txt");
    convergence_table.write_text(output_1, dealii::TableHandler::TextOutputFormat::org_mode_table);

    convergence_table.add_column_to_supercolumn("cycle", "n cells");
    convergence_table.add_column_to_supercolumn("cells", "n cells");
    std::vector<std::string> new_order;
    new_order.emplace_back("n cells");
    new_order.emplace_back("L2");
    convergence_table.set_column_order(new_order);
    if (true && parameter_set->refine_time)
        convergence_table.evaluate_convergence_rates(
            "L2", dealii::ConvergenceTable::reduction_rate_log2);

    // compute convergence rates from 3 consecutive errors
    for (unsigned int i = 0; i < L2_error_vals.size(); ++i) {
        if (i < 2)
        convergence_table.add_value("L2...", std::string("-"));
        else {
        double p0 = L2_error_vals[i - 2];
        double p1 = L2_error_vals[i - 1];
        double p2 = L2_error_vals[i];
        convergence_table.add_value(
            "L2...", std::log(std::fabs((p0 - p1) / (p1 - p2))) / std::log(2.));
        }
    }

    std::ofstream output(output_dir + "convergence.txt");
    convergence_table.write_text(output, dealii::TableHandler::TextOutputFormat::org_mode_table);
    std::cout << std::endl;
    convergence_table.write_text(std::cout, dealii::TableHandler::TextOutputFormat::org_mode_table);
}

template <int dim>
void ElasticWave_cGcG<dim>::run() {
    init_functions();
    init_grid();
    
    refinement_table.print_header();
    refinement_table.add_value("total", 4);


    for (unsigned int cycle = 0; cycle < parameter_set->max_refinements_time; ++cycle) {
        refinement_table.add_value("refinement cycle", cycle);
        refinement_table.print_row();
        
        const std::string output_dir = "../output/cG-cG/example-" + std::to_string(parameter_set->example) + "/"
                                    + "final_condition=" + std::to_string(parameter_set->junker) + "/" 
                                    + "cylce=" + std::to_string(cycle) + "/";
        
        mkdir(output_dir.c_str(), 0777);

        L2_error = 0.;      

        dealii::Timer timer;

        setup_dofs();
        refinement_table.add_value("current slab", 1);
        refinement_table.print_row();
        assemble_system();
        refinement_table.add_value("current slab", 2);
        refinement_table.print_row();
        solve();
        refinement_table.add_value("current slab", 3);
        refinement_table.print_row();
        output_results(cycle);
        refinement_table.add_value("current slab", 4);
        refinement_table.print_row();
        timer.stop();

        // if exact solution exists error is processed
        if (function.exact_solution != nullptr)
            process_solution(cycle);
        

        refinement_table.print_row();
        
        refinement_table.add_value("time", timer.wall_time());
        refinement_table.print_row();
        
        
        if (dim == 1)
            plot_solution(cycle);
        grid.refine_global(parameter_set->refine_space, parameter_set->refine_time);


    }
    
    if (function.exact_solution != nullptr)
        print_convergence_table();

    
    
}

#include "ElasticWave_cGcG.inst.in"
