//*---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the ideal.II authors
//
// This file is part of the ideal.II library.
//
// The ideal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of ideal.II.
//
//*---------------------------------------------------------------------




// PROJECT includes
#include <wave/ElasticWave_dG.tpl.hh>

#include <wave/InitialValue/InitialValue_Selector.tpl.hh>
#include <wave/FinalValue/FinalValue_Selector.tpl.hh>
#include <wave/Force/Force_Selector.tpl.hh>
#include <wave/Boundary/BoundaryValue_Selector.tpl.hh>
#include <wave/ExactSolution/ExactSolution_Selector.tpl.hh>
#include <wave/ElasticWave_dG.tpl.hh>


// iDEAL.II includes
#include <ideal.II/numerics/vector_tools.hh>


// DEAL.II includes



// C++ includes




template <int dim>
ElasticWave_dG<dim>::ElasticWave_dG(unsigned int s, unsigned int r)
    :   triangulation(1),
        dof_handler(&triangulation),
        fe(std::make_shared< dealii::FESystem<dim> >(/*u*/ dealii::FE_Q<dim>(s), dim,
                                                     /*v*/ dealii::FE_Q<dim>(s), dim), 
           r),
        slab(0)
{}

template <int dim>
inline dealii::SymmetricTensor<2, dim> ElasticWave_dG<dim>::get_strain(const dealii::Tensor<2,dim> &grad) {
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
void ElasticWave_dG<dim>::set_input_parameters(
    std::shared_ptr< dealii::ParameterHandler > parameter_handler) {
    
    parameter_set = std::make_shared< wave::ParameterSet > (
		parameter_handler); 
}


template <int dim>
void ElasticWave_dG<dim>::init_functions() {


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
void ElasticWave_dG<dim>::init_grid() {

    // construct a shared pointer to a spatial triangulation
    auto space_tria = std::make_shared<dealii::Triangulation<dim>>();

    // generator to generate spatial triangulation
    wave::grid::TriaGenerator<dim> tria_gen;
    tria_gen.generate(space_tria, parameter_set->example);

    // fill the list with parameter_set->M ""slab::Triangulation"" objects sharing the same spatial triangulation
    triangulation.generate(space_tria, parameter_set->M, parameter_set->t0, parameter_set->T);
    // refine the grids on each slab in (space, time)
    triangulation.refine_global(parameter_set->global_refinement_space, 0);
    // generate a slab::DoFHandler for each slab::Triangulation
    dof_handler.generate();
    
}

template <int dim>
void ElasticWave_dG<dim>::time_steps(const unsigned int refinement_cycle) {

    // This collection simplifies the time step increments
    idealii::TimeIteratorCollection<dim> it_coll = idealii::TimeIteratorCollection<dim>();

    // Fill the solution list of vectors with n_slabs dealii::Vector<double>
    solution.reinit(triangulation.M());

    // Get first iterators to the first slab
    slab_it.tria = triangulation.begin();
    slab_it.dof = dof_handler.begin();
    slab_it.solution = solution.begin();

    // Register to the TimeIteratorCollector
    it_coll.add_iterator(&slab_it.tria, &triangulation);
    it_coll.add_iterator(&slab_it.dof, &dof_handler);
    it_coll.add_iterator(&slab_it.solution, &solution);

    refinement_table.add_value("total", triangulation.M());

    slab = 0;
    dealii::Timer timer;

    // actual Time Steps using ""increment""
    for (; !it_coll.at_end(); it_coll.increment()) {
                
        // typical FEM loop without ""inite_grid""
        setup_dofs_on_slab();
        assemble_system_on_slab();
        solve_on_slab();
        output_results(refinement_cycle);

        // if exact solution exists error is processed
        if (function.exact_solution != nullptr)
            process_solution(refinement_cycle);
        
        // extract the subvector of the final DoF of this slab for initial value of the next slab
        idealii::slab::VectorTools::extract_subvector_at_time_dof(*slab_it.solution,
                                                                  slab_initial_value,
                                                                  slab_it.tria->temporal()->n_global_active_cells() - 1);


        refinement_table.add_value("current slab", slab + 1);
        refinement_table.add_value("refinement cycle", refinement_cycle);
        const unsigned int n_dofs_space = slab_it.dof->n_dofs_space();
        const unsigned int n_dofs_time  = slab_it.dof->n_dofs_time() * triangulation.M();
        refinement_table.add_value("dofs (space)", n_dofs_space);
        refinement_table.add_value("dofs (time)", n_dofs_time);

        refinement_table.print_row();
        // increase slab index
        slab++;
        
    }

    timer.stop();
    refinement_table.add_value("time", timer.wall_time());
    refinement_table.print_row();
    
}


template <int dim>
void ElasticWave_dG<dim>::setup_dofs_on_slab() {

    // Distribute spatial and temporal dofs
    slab_it.dof->distribute_dofs(fe);

    std::vector<unsigned int> block_component(dim + dim, 0);
    std::fill(block_component.begin() + dim, block_component.end(), 1);


    // Renumbering spartial DoFs into displacement and velocity
    dealii::DoFRenumbering::component_wise(*slab_it.dof->spatial(), block_component);
       
    // two blocks: displacement and velocity
    const std::vector<dealii::types::global_dof_index> dofs_per_block = 
          dealii::DoFTools::count_dofs_per_fe_block(*slab_it.dof->spatial(), block_component);

    n_space_u = dofs_per_block[0];
    n_space_v = dofs_per_block[1];


    // on the first slab the given initial value is used
    if (slab == 0) {
        // computing Initial Values
        slab_initial_value.reinit(slab_it.dof->n_dofs_space());
        dealii::VectorTools::interpolate(*(slab_it.dof->spatial()), *function.inital_value_function, slab_initial_value, dealii::ComponentMask());
        
    }

    // set homogenouse Dirichlet boundary constraints
    slab_constraints = std::make_shared< dealii::AffineConstraints<double> > ();
    idealii::slab::VectorTools::interpolate_boundary_values(*slab_it.dof, 0,*function.boundary_function, slab_constraints);
    slab_constraints->close();


    // Construct the space-time sparsity pattern for this slab
    dealii::DynamicSparsityPattern dsp(slab_it.dof->n_dofs_spacetime());
    idealii::slab::DoFTools::make_upwind_sparsity_pattern(*slab_it.dof, dsp);
    slab_sparsity_pattern.copy_from(dsp);

    // print sparsity_pattern of the space-time slab
    // // std::ofstream out_sparsity("../sparsity_pattern.svg");
    // // slab_sparsity_pattern.print_svg(out_sparsity);

    // reinit the linear system
    slab_system_matrix.reinit(slab_sparsity_pattern);
    slab_it.solution->reinit(slab_it.dof->n_dofs_spacetime());
    slab_system_rhs.reinit(slab_it.dof->n_dofs_spacetime());

}


template <int dim>
void ElasticWave_dG<dim>::assemble_system_on_slab() {



    // reset slab system matrix and right hand side
    // slab_system_matrix = 0;
    // slab_system_rhs = 0;

    // space - time quadratur and FEValues
    idealii::spacetime::QGauss<dim> quad(fe.spatial()->degree + 3, 
                                         fe.temporal()->degree + 3);
    idealii::spacetime::FEValues<dim> fe_values(fe, quad,
                                                dealii::update_values | dealii::update_gradients | 
                                                dealii::update_quadrature_points | dealii::update_JxW_values);

    // To account the jump values we use the ""FEJumpValues"" class (similar to the FEFaceValues)
    idealii::spacetime::FEJumpValues<dim> fe_jump_values(fe, quad,
                                                         dealii::update_values | dealii::update_gradients | 
                                                         dealii::update_quadrature_points | dealii::update_JxW_values);

    const unsigned int dofs_per_spacetime_cell = fe.dofs_per_cell;

    // Number of temporal elements and index of the current elelemt for offset calculations
    unsigned int N = slab_it.tria->temporal()->n_global_active_cells();
    unsigned int n;

    // space-time local matrix and rhs containing all temporal information but only one spatial element
    dealii::FullMatrix<double> local_matrix(N * dofs_per_spacetime_cell, 
                                            N * dofs_per_spacetime_cell);
    dealii::Vector<double> local_rhs(N * dofs_per_spacetime_cell);
    std::vector<dealii::types::global_dof_index> local_dof_indices(N * 
                                                                   dofs_per_spacetime_cell);

    // number of space-time and space quadratur points
    unsigned int n_quad_spacetime = fe_values.n_quadrature_points;
    unsigned int n_quad_space = quad.spatial()->size();

    double lame_coefficient_lambda = (2 * parameter_set->poisson_ratio_nu * parameter_set->lame_coefficient_mu) / (1.0 - 2 * parameter_set->poisson_ratio_nu);
    
    dealii::SymmetricTensor<2, dim> identity = dealii::unit_symmetric_tensor<dim, double>();
  
    // First iterate over active spatial elements
    for (const auto &space_cell : slab_it.dof->spatial()->active_cell_iterators()) {

        // reset local contributions
        local_matrix = 0;
        local_rhs = 0;

        // recalculate local information for the current spatial element
        fe_values.reinit_space(space_cell);
        fe_jump_values.reinit_space(space_cell);

        // get local contribution of the slab_initial_vlaue vector
        // displacement
        std::vector<dealii::Tensor<1, dim>> initial_value_u(fe_values.spatial()->n_quadrature_points);
        fe_values.spatial()->operator[](displacement).get_function_values(slab_initial_value, initial_value_u); 
       
        // velocity
        std::vector<dealii::Tensor<1, dim>> initial_value_v(fe_values.spatial()->n_quadrature_points);
        fe_values.spatial()->operator[](velocity).get_function_values(slab_initial_value, initial_value_v);     


        // iterate over all active temporal elements
        for (const auto &time_cell : slab_it.dof->temporal()->active_cell_iterators()) {
            n = time_cell->index();
            
            // recalculate the local information of the current temporla element
            fe_values.reinit_time(time_cell);
            fe_jump_values.reinit_time(time_cell);

            // Get local space-time dof-indices for the current tensor product cell
            fe_values.get_local_dof_indices(local_dof_indices);
                        
            // Iterate over all quadratur points except for jump terms
            for (unsigned int q = 0; q < n_quad_spacetime; ++q) {

                // set the time of the right hand side of the current quadratur point and get its location
                const double t_qq = fe_values.time_quadrature_point(q);
                function.rhs_function->set_time(t_qq);
                const auto &x_q = fe_values.space_quadrature_point(q);

                // Iterate over all space-time dofs of the current element
                for (unsigned int i = 0; i < dofs_per_spacetime_cell; ++i) {

                    local_rhs(i + n * dofs_per_spacetime_cell) += 
                        (fe_values.vector_value(displacement, i, q) *                       // φ^u_{i_time, i_space}(t_qq, x_q)
                         function.rhs_function->value(x_q) *                                // f(t_qq, x_q)
                         fe_values.JxW(q));                                                 // dx dt


                    for (unsigned int j = 0; j < dofs_per_spacetime_cell; ++j) {

                        local_matrix(i + n * dofs_per_spacetime_cell,
                                     j + n * dofs_per_spacetime_cell) += 
                            
                             (parameter_set->rho *                                          // ρ
                              fe_values.vector_value(displacement, i, q) *                  // φ^u_{j_time, j_space}(t_qq, x_q)
                              fe_values.vector_dt(velocity, j, q)                           // ∂_t φ^v_{i_time, i_space}(t_qq, x_q)


                            + 2 * parameter_set->lame_coefficient_mu *                          // 2 * μ
                              dealii::scalar_product(
                              fe_values.vector_space_grad(displacement, i, q),                  // ∇_x φ^u_{i_time, i_space}(t_qq, x_q)
                              get_strain(fe_values.vector_space_grad(displacement, j, q)))      // 1/2 * (∇_x φ^u_{j_time, j_space}(t_qq, x_q) + ∇_x^T φ^u_{j_time, j_space}(t_qq, x_q))
                                 
                                
                            + lame_coefficient_lambda *                                          // λ
                              dealii::scalar_product(
                              fe_values.vector_space_grad(displacement, i, q),                                          // ∇_x φ^u_{i_time, i_space}(t_qq, x_q)
                              dealii::trace(get_strain(fe_values.vector_space_grad(displacement, j, q))) * identity)    // tr(1/2 * (∇_x φ^u_{j_time, j_space}(t_qq, x_q) + ∇_x^T φ^u_{j_time, j_space}(t_qq, x_q)))
                                                           
                                                          
                            + fe_values.vector_value(velocity, i, q) *                      // φ^v_{j_time, j_space}(t_qq, x_q)
                              fe_values.vector_dt(displacement, j, q)                       // ∂_t φ^u_{i_time, i_space}(t_qq, x_q)


                            - fe_values.vector_value(velocity, i, q) *                      // φ^v_{i_time, i_space}(t_qq, x_q)
                              fe_values.vector_value(velocity, j, q))                       // φ^v_{j_time, j_space}(t_qq, x_q)

                            * fe_values.JxW(q);                                             // dx dt
                    } // dof j
                } // dof i
                    
            } // q



            // Jump terms if dG is used for the time
           
                // Jump terms only have a spatial quadratur loop
                for (unsigned int q = 0; q < n_quad_space; ++q) {
                    for (unsigned int i = 0; i < dofs_per_spacetime_cell; ++i) {
                        for (unsigned int j = 0; j < dofs_per_spacetime_cell; ++j) {
                            local_matrix(i + n * dofs_per_spacetime_cell,
                                        j + n * dofs_per_spacetime_cell) += 
                                
                                  (fe_jump_values.vector_value_plus(velocity, i, q) *           // φ^{v,+}_{j_time, j_space}(t_qq^+, x_q)
                                   fe_jump_values.vector_value_plus(displacement, j, q)         // φ^{u,+}_{i_time, i_space}(t_qq^+, x_q)

                                + fe_jump_values.vector_value_plus(displacement, i, q) *        // φ^{v,+}_{j_time, j_space}(t_qq^+, x_q)
                                  fe_jump_values.vector_value_plus(velocity, j, q))             // φ^{u,+}_{i_time, i_space}(t_qq^+, x_q)
                                
                                * fe_jump_values.JxW(q);                                        // dx 

                            
                            if (n > 0) {
                                local_matrix(i + n * dofs_per_spacetime_cell,
                                        j + (n - 1) * dofs_per_spacetime_cell) -= 
                                
                                  (fe_jump_values.vector_value_plus(velocity, i, q) *           // φ^{v,+}_{j_time, j_space}(t_qq^+, x_q)
                                   fe_jump_values.vector_value_minus(displacement, j, q)        // φ^{u,-}_{i_time, i_space}(t_qq^-, x_q)

                                + fe_jump_values.vector_value_plus(displacement, i, q) *        // φ^{u,+}_{j_time, j_space}(t_qq^+, x_q)
                                  fe_jump_values.vector_value_minus(velocity, j, q))            // φ^{v,-}_{i_time, i_space}(t_qq^-, x_q)
                                
                                * fe_jump_values.JxW(q);                                        // dx 
                            }
                        } // dof j

                        if (n == 0) {
                            local_rhs(i) += 
                                  (fe_jump_values.vector_value_plus(velocity,i, q) *             // φ^{v,+}_{i_time, i_space}(0^+, x_q)
                                   initial_value_u[q]                                            // u_0(x_q)
                                  

                                + fe_jump_values.vector_value_plus(displacement,i, q) *         // φ^{u,+}_{i_time, i_space}(0^+, x_q)
                                  initial_value_v[q])                                           // v_0(x_q)
                                
                                * fe_jump_values.JxW(q);                                        // dx
                        }


                    } // dof i

                } // q


        } // cell time

        
        slab_constraints->distribute_local_to_global(local_matrix, local_rhs, local_dof_indices, slab_system_matrix, slab_system_rhs);


    } // cell space



}

template <int dim>
void ElasticWave_dG<dim>::solve_on_slab() {


    dealii::SparseDirectUMFPACK solver;
    solver.factorize(slab_system_matrix);
    solver.vmult(*slab_it.solution, slab_system_rhs);


    slab_constraints->distribute(*slab_it.solution);

}

template <int dim>
void ElasticWave_dG<dim>::output_results(const unsigned int refinement_cycle) {
    
    const std::string output_dir = "../output/dG/example-" + std::to_string(parameter_set->example) + "/"
                                    + "final_condition=" + std::to_string(parameter_set->junker) + "/" 
                                    + "cylce=" + std::to_string(refinement_cycle) + "/";

    std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> 
            data_component_interpretation(2*dim, dealii::DataComponentInterpretation::component_is_part_of_vector);
    std::vector<std::string> solution_names(dim, "displacement");
    std::vector<std::string> solution_names_v(dim, "velocity");
    solution_names.insert(solution_names.end(), solution_names_v.begin(), solution_names_v.end());

    auto n_dofs = slab_it.dof->n_dofs_time();

    // time
    std::vector<dealii::Point<1>> time_support_points(slab_it.dof->temporal()->n_dofs());
    dealii::DoFTools::map_dofs_to_support_points(
        dealii::MappingQ1<1, 1>(), *slab_it.dof->temporal(), time_support_points);

    
    for (unsigned int i = 0; i < n_dofs; i++) {

        dealii::Point<1> point = time_support_points[i];
        
        dealii::DataOut<dim> data_out;
        data_out.attach_dof_handler(*slab_it.dof->spatial());
        dealii::Vector<double> local_solution;
        local_solution.reinit(slab_it.dof->n_dofs_space());

        idealii::slab::VectorTools::extract_subvector_at_time_dof(*slab_it.solution, local_solution, i);
     
        data_out.add_data_vector(local_solution, solution_names, dealii::DataOut<dim>::type_dof_data, data_component_interpretation);
        data_out.build_patches();
        data_out.set_flags(dealii::DataOutBase::VtkFlags(point[0]));
        
        std::ostringstream filename;
        filename << "solution(" << fe.temporal()->degree << ")_t_" << slab * n_dofs + i << ".vtk";

        std::ofstream output(output_dir + filename.str());
        data_out.write_vtk(output);

        output.close();

    }    
    
}


template <int dim>
void ElasticWave_dG<dim>::process_solution(unsigned int cycle) {


    idealii::spacetime::QGauss<dim> quad(fe.spatial()->degree + 2, fe.temporal()->degree + 2);
    L2_error += idealii::slab::VectorTools::calculate_L2L2_squared_error_on_slab<dim>(*slab_it.dof, *slab_it.solution, *function.exact_solution, quad);


    if (slab == (triangulation.M() - 1)) {
        L2_error = std::sqrt(L2_error);
        L2_error_vals.push_back(L2_error);

        const unsigned int n_dofs_space = slab_it.dof->n_dofs_space();
        const unsigned int n_dofs_time  = slab_it.dof->n_dofs_time() * triangulation.M();
        const unsigned int n_dofs = n_dofs_space * n_dofs_time;

        convergence_table.add_value("cycle", cycle);
        convergence_table.add_value("cells", slab_it.tria->spatial()->n_active_cells() * triangulation.M());
        convergence_table.add_value("dofs", n_dofs);
        convergence_table.add_value("dofs(space)", n_dofs_space);
        convergence_table.add_value("dofs(time)", n_dofs_time);
        convergence_table.add_value("L2", L2_error);
    }
    
}





template <int dim>
void ElasticWave_dG<dim>::plot_solution(const unsigned int refinement_cycle) {

    const std::string output_dir = "../output/dG/example-" + std::to_string(parameter_set->example) + "/"
                                    + "final_condition=" + std::to_string(parameter_set->junker) + "/" 
                                    + "cylce=" + std::to_string(refinement_cycle) + "/";

    std::string py_cmd2 = "python3 ../python/plot_solution.py " + output_dir + " " + std::to_string(parameter_set->example);
    system(py_cmd2.c_str());
}

template <int dim>
void ElasticWave_dG<dim>::print_convergence_table() {

    const std::string output_dir = "../output/dG/example-" + std::to_string(parameter_set->example) + "/"
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
    // if (true && parameter_set->refine_time)
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
void ElasticWave_dG<dim>::run() {
    init_functions();
    init_grid();
    
    
    refinement_table.print_header();


    for (unsigned int cycle = 0; cycle < parameter_set->max_refinements_time; ++cycle) {
        
        
        const std::string output_dir = "../output/dG/example-" + std::to_string(parameter_set->example) + "/"
                                    + "final_condition=" + std::to_string(parameter_set->junker) + "/" 
                                    + "cylce=" + std::to_string(cycle) + "/";
        
        mkdir(output_dir.c_str(), 0777);

        L2_error = 0.;

        time_steps(cycle);
        
        
        if (dim == 1)
            plot_solution(cycle);
        triangulation.refine_global(parameter_set->refine_space, parameter_set->refine_time);
        
        // reinit a slab::DoFHandler for each slab::Triangulation
        dof_handler.reinit();

    }
    
    if (function.exact_solution != nullptr)
        print_convergence_table();
    
    
}

// Implementierung für die Dimension 1
#include "ElasticWave_dG.inst.in"
