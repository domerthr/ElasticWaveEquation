#ifndef __vector_tools_hh
#define __vector_tools_hh

// PROJECT includes
#include <fest/base/quadrature_lib.tpl.hh>
#include <fest/dofs/SpaceTime_dof_handler.tpl.hh>
#include <fest/dofs/Slab_dof_handler.tpl.hh>

// iDEAL.II includes
#include <ideal.II/base/quadrature_lib.hh>


// DEAL.II includes
#include <deal.II/base/config.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>


// C++ includes
#include <memory>



namespace fest {
namespace cgcg {
namespace VectorTools {



template <int dim, typename Number>
void interpolate_boundary_values(fest::cgcg::DoFHandler<dim> &dof_handler,
                                 const dealii::types::boundary_id boundary_component,
                                 dealii::Function<dim, Number> &boundary_function,
                                 std::shared_ptr<dealii::AffineConstraints<Number>> spacetime_constraints,
                                 const dealii::ComponentMask &component_mask = dealii::ComponentMask());

template <int dim>
double calculate_L2L2_squared_error(fest::cgcg::DoFHandler<dim> &dof_handler,
                                    dealii::Vector<double> &space_time_vector,
                                    dealii::Function<dim, double> &exact_solution,
                                    fest::Quadrature<dim> &quad);


void extract_subvector_at_time_dof(const dealii::Vector<double> &spacetime_vector,
                                   dealii::Vector<double> &space_vector,
                                   unsigned int dof_index);



} // namespace VectorTools
} // namespace cgcg
namespace cg {
namespace VectorTools {

template <int dim, typename Number>
void interpolate_boundary_values(fest::cg::slab::DoFHandler<dim> &dof_handler,
                                 const dealii::types::boundary_id boundary_component,
                                 dealii::Function<dim, Number> &boundary_function,
                                 std::shared_ptr<dealii::AffineConstraints<Number>> spacetime_constraints,
                                 const dealii::ComponentMask &component_mask = dealii::ComponentMask());

template <int dim>
double calculate_L2L2_squared_error_on_slab(fest::cg::slab::DoFHandler<dim> &dof_handler,
                                            dealii::Vector<double> &space_time_vector,
                                            dealii::Function<dim, double> &exact_solution,
                                            idealii::spacetime::Quadrature<dim> &quad);


void extract_subvector_at_time_dof(const dealii::Vector<double> &spacetime_vector,
                                   dealii::Vector<double> &space_vector,
                                   unsigned int dof_index);



} // namespace VectorTools
} // namespace cg
} // namespace fest

#endif
