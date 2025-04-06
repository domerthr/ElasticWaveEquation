// PROJECT includes
#include <fest/numerics/vector_tools.hh>


// DEAL.II includes
#include <deal.II/dofs/dof_tools.h>



// C++ includes


namespace fest {
namespace cgcg {
namespace VectorTools {

template <int dim, typename Number>
void interpolate_boundary_values(fest::cgcg::DoFHandler<dim> &dof_handler,
                                 const dealii::types::boundary_id boundary_component,
                                 dealii::Function<dim, Number> &boundary_function,
                                 std::shared_ptr< dealii::AffineConstraints<Number> > spacetime_constraints,
                                 const dealii::ComponentMask &component_mask) {

    auto space_constraints = std::make_shared<dealii::AffineConstraints<Number>>();

    dealii::Quadrature<1> quad_time(dof_handler.temporal()->get_fe(0).get_unit_support_points());
    dealii::FEValues<1> fev(dof_handler.temporal()->get_fe(0),
                            quad_time,
                            dealii::update_quadrature_points);

    unsigned int n_space_dofs = dof_handler.n_dofs_space();
    unsigned int time_dof     = 0;
    std::vector<dealii::types::global_dof_index> local_indices(dof_handler.dofs_per_cell_time());
    // loop over time cells instead
    for (const auto &cell_time : dof_handler.temporal()->active_cell_iterators()) {

        fev.reinit(cell_time);

        cell_time->get_dof_indices(local_indices);
        for (unsigned int q = 0; q < quad_time.size(); q++) {

            space_constraints->clear();
            time_dof = local_indices[q];
            boundary_function.set_time(fev.quadrature_point(q)[0]);
            dealii::VectorTools::interpolate_boundary_values(*dof_handler.spatial(), boundary_component,
                                                             boundary_function, *space_constraints, component_mask);

            for (unsigned int i = 0; i < n_space_dofs; i++) {
                // check if this is a constrained dof
                if (space_constraints->is_constrained(i)) {

                    const std::vector<std::pair<dealii::types::global_dof_index, double>> *entries = space_constraints->get_constraint_entries(i);

                    spacetime_constraints->add_line(i + time_dof * n_space_dofs);

                    // non Dirichlet constraint
                    if (entries->size() > 0) {
                        for (auto entry : *entries) {
                            spacetime_constraints->add_entry(i + time_dof * n_space_dofs, entry.first + time_dof * n_space_dofs, entry.second);
                        }
                    }
                    else{

                        spacetime_constraints->set_inhomogeneity(i + time_dof * n_space_dofs, space_constraints->get_inhomogeneity(i));
                    }
                }
            }
        }
    }
}



template <int dim, typename Number>
void interpolate_boundary_values(dealii::IndexSet space_relevant_dofs,
                                 fest::cgcg::DoFHandler<dim> &dof_handler,
                                 const dealii::types::boundary_id boundary_component,
                                 dealii::Function<dim, Number> &boundary_function,
                                 std::shared_ptr<dealii::AffineConstraints<Number>> spacetime_constraints,
                                 const dealii::ComponentMask &component_mask = dealii::ComponentMask())
  {
    auto space_constraints =
      std::make_shared<dealii::AffineConstraints<Number>>();
    space_constraints->reinit(space_relevant_dofs);
    dealii::Quadrature<1> quad_time(
      dof_handler.temporal()->get_fe(0).get_unit_support_points());
    dealii::FEValues<1> fev(dof_handler.temporal()->get_fe(0),
                            quad_time,
                            dealii::update_quadrature_points);

    unsigned int n_space_dofs = dof_handler.n_dofs_space();
    unsigned int time_dof     = 0;
    std::vector<dealii::types::global_dof_index> local_indices(
      dof_handler.dofs_per_cell_time());
    // loop over time cells instead
    for (const auto &cell_time : dof_handler.temporal()->active_cell_iterators())
    {
        fev.reinit(cell_time);

        cell_time->get_dof_indices(local_indices);
        for (unsigned int q = 0; q < quad_time.size(); q++) {
            space_constraints->clear();
            time_dof = local_indices[q];
            boundary_function.set_time(fev.quadrature_point(q)[0]);
            dealii::VectorTools::interpolate_boundary_values(
              *dof_handler.spatial(),
              boundary_component,
              boundary_function,
              *space_constraints,
              component_mask);

            for (auto id = space_relevant_dofs.begin();
                 id != space_relevant_dofs.end();
                 id++) {
                // check if this is a constrained dof
                if (space_constraints->is_constrained(*id)) {
                    const std::vector<std::pair<dealii::types::global_dof_index,
                                                double>> *entries =
                    space_constraints->get_constraint_entries(*id);
                    spacetime_constraints->add_line(*id +
                                                    time_dof * n_space_dofs);
                    // non Dirichlet constraint
                    if (entries->size() > 0) {
                        for (auto entry : *entries) {
                            std::cout << entry.first << "," << entry.second
                                      << std::endl;
                            spacetime_constraints->add_entry(
                              *id + time_dof * n_space_dofs,
                              entry.first + time_dof * n_space_dofs,
                              entry.second);
                        }
                    }
                    else {
                        spacetime_constraints->set_inhomogeneity(
                          *id + time_dof * n_space_dofs,
                          space_constraints->get_inhomogeneity(*id));
                    }
                }
            }
        }
    }
}




template <int dim>
double calculate_L2L2_squared_error(fest::cgcg::DoFHandler<dim> &dof_handler,
                                    dealii::Vector<double> &space_time_vector,
                                    dealii::Function<dim, double> &exact_solution,
                                    fest::Quadrature<dim> &quad) {

    double L2_error_squared = 0.;
    dealii::FEValues<1> time_fe_values(dof_handler.temporal()->get_fe(),
                                       *quad.temporal(),
                                       dealii::update_values | dealii::update_JxW_values | dealii::update_quadrature_points);
    
    unsigned int n_dofs_space = dof_handler.n_dofs_space();
    
    dealii::Vector<double> space_solution;
    dealii::Vector<double> difference_per_cell;

    space_solution.reinit(n_dofs_space);
    difference_per_cell.reinit(n_dofs_space);

    std::vector< dealii::Point<1,double> > q_points;
    double t = 0;
    unsigned int offset = 0;

    for (const auto &time_cell : dof_handler.temporal()->active_cell_iterators()) {

        offset = time_cell->index() *  n_dofs_space;

        time_fe_values.reinit(time_cell);

        for (unsigned int q = 0; q < time_fe_values.n_quadrature_points; q++) {
            // time point
            t = time_fe_values.quadrature_point(q)[0];
            exact_solution.set_time(t);

            space_solution = 0;

            for (unsigned int ii = 0; ii < dof_handler.dofs_per_cell_time(); ii++) {
                double factor = time_fe_values.shape_value(ii, q);
                for (unsigned int i = 0; i < n_dofs_space; i++) {

                    space_solution[i] += space_time_vector[i + ii * n_dofs_space + offset] * factor;
                }

            }

        
            difference_per_cell = 0;
            dealii::VectorTools::integrate_difference(*dof_handler.spatial(), space_solution, exact_solution, difference_per_cell,
                                                      *quad.spatial(), dealii::VectorTools::L2_norm);
            

            // add local contributions to global L2 error
            L2_error_squared += difference_per_cell.norm_sqr() * time_fe_values.JxW(q);
        
        }

    }
    

    return L2_error_squared;

}



void extract_subvector_at_time_dof(const dealii::Vector<double> &spacetime_vector,
                                   dealii::Vector<double> &space_vector,
                                   unsigned int dof_index) {

    unsigned int n_dofs_space = space_vector.size();
    for (unsigned int i = 0; i < n_dofs_space; i++) {
        space_vector[i] = spacetime_vector[i + dof_index * n_dofs_space];
    }
}


} // namespace VectorTools
} // namespace cgcg


namespace cg {
namespace VectorTools {

template <int dim, typename Number>
void
interpolate_boundary_values(fest::cg::slab::DoFHandler<dim>                    &dof_handler,
                            const dealii::types::boundary_id                   boundary_component,
                            dealii::Function<dim, Number>                      &boundary_function,
                            std::shared_ptr<dealii::AffineConstraints<Number>> spacetime_constraints,
                            const dealii::ComponentMask                        &component_mask)
{
    auto space_constraints =
    std::make_shared<dealii::AffineConstraints<Number>>();

    dealii::Quadrature<1> quad_time(
        dof_handler.temporal_trial()->get_fe(0).get_unit_support_points());
    dealii::FEValues<1> fev(dof_handler.temporal_trial()->get_fe(0),
                            quad_time,
                            dealii::update_quadrature_points);

    unsigned int n_space_dofs = dof_handler.n_dofs_space();
    unsigned int time_dof     = 0;
    std::vector<dealii::types::global_dof_index> local_indices(
        dof_handler.dofs_per_cell_time_trial());
    // loop over time cells instead
    for (const auto &cell_time : dof_handler.temporal_trial()->active_cell_iterators()) {
        fev.reinit(cell_time);

        cell_time->get_dof_indices(local_indices);
        for (unsigned int q = 0; q < quad_time.size(); q++) {

            space_constraints->clear();
            time_dof = local_indices[q];
            boundary_function.set_time(fev.quadrature_point(q)[0]);
            dealii::VectorTools::interpolate_boundary_values(
                *dof_handler.spatial(),
                boundary_component,
                boundary_function,
                *space_constraints,
                component_mask);

            for (unsigned int i = 0; i < n_space_dofs; i++) {

                // check if this is a constrained dof
                if (space_constraints->is_constrained(i))
                {
                    const std::vector<std::pair<dealii::types::global_dof_index,
                                                double>> *entries =
                        space_constraints->get_constraint_entries(i);

                    spacetime_constraints->add_line(i +
                                                    time_dof * n_space_dofs);
                    // non Dirichlet constraint
                    if (entries->size() > 0)
                    {
                        for (auto entry : *entries)
                        {
                            spacetime_constraints->add_entry(
                                i + time_dof * n_space_dofs,
                                entry.first + time_dof * n_space_dofs,
                                entry.second);
                        }
                    }
                    else
                    {
                        spacetime_constraints->set_inhomogeneity(
                            i + time_dof * n_space_dofs,
                            space_constraints->get_inhomogeneity(i));
                    }
                }
            }
        }
    }
}

template <int dim>
double
calculate_L2L2_squared_error_on_slab(fest::cg::slab::DoFHandler<dim>         &dof_handler,
                                     dealii::Vector<double>        &spacetime_vector,
                                     dealii::Function<dim, double> &exact_solution,
                                     idealii::spacetime::Quadrature<dim>    &quad)
{
    double              slab_norm_sqr = 0.;
    dealii::FEValues<1> fev_time(dof_handler.temporal_trial()->get_fe(),
                                 *quad.temporal(),
                                 dealii::update_values |
                                 dealii::update_JxW_values |
                                 dealii::update_quadrature_points);


    dealii::Vector<double> space_vec;
    dealii::Vector<double> difference_per_cell;
    space_vec.reinit(dof_handler.n_dofs_space());
    difference_per_cell.reinit(dof_handler.n_dofs_space());
    std::vector<dealii::Point<1, double>> q_points;
    double                                t      = 0;
    unsigned int                          offset = 0;
    unsigned int n_dofs_space                    = dof_handler.n_dofs_space();
    for (auto cell_time : dof_handler.temporal_trial()->active_cell_iterators())
    {
        offset =
            cell_time->index() * dof_handler.dofs_per_cell_time_trial() * n_dofs_space;

        fev_time.reinit(cell_time);
        for (unsigned int q = 0; q < fev_time.n_quadrature_points; q++)
        {
            t = fev_time.quadrature_point(q)[0];
            exact_solution.set_time(t);
            // calculate spatial vector at t
            space_vec = 0;
            for (unsigned int ii = 0; ii < dof_handler.dofs_per_cell_time_trial();
                    ii++)
            {
                double factor = fev_time.shape_value(ii, q);
                for (unsigned int i = 0; i < n_dofs_space; i++)
                {
                    space_vec[i] +=
                        spacetime_vector[i + ii * n_dofs_space + offset] * factor;
                }
            }
            // calculate L2 norm at current temporal QP
            difference_per_cell = 0;
            dealii::VectorTools::integrate_difference(
                *dof_handler.spatial(),
                space_vec,
                exact_solution,
                difference_per_cell,
                *quad.spatial(),
                dealii::VectorTools::L2_norm);

            slab_norm_sqr += difference_per_cell.norm_sqr() * fev_time.JxW(q);
        }
    }
    return slab_norm_sqr;
}


void extract_subvector_at_time_dof(const dealii::Vector<double> &spacetime_vector,
                                   dealii::Vector<double> &space_vector,
                                   unsigned int dof_index) {

    unsigned int n_dofs_space = space_vector.size();
    for (unsigned int i = 0; i < n_dofs_space; i++) {
        space_vector[i] = spacetime_vector[i + dof_index * n_dofs_space];
    }
}

} // namespace VectorTools
} // namespace cg


} // namespace fest


#include "vector_tools.inst.in"