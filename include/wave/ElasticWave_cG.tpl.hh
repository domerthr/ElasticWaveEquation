#ifndef __ElasticWave_cG_tpl_hh
#define __ElasticWave_cG_tpl_hh

// PROJECT includes
#include <fest/grid/SpaceTimeTria.tpl.hh>

#include <fest/dofs/SpaceTime_dof_handler.tpl.hh>
#include <fest/dofs/Slab_dof_handler.tpl.hh>
#include <fest/dofs/SpaceTime_dof_tools.hh>

#include <fest/fe/fe.tpl.hh>
#include <fest/fe/SpaceTime_fe_values.tpl.hh>

#include <fest/base/refinement_table.hh>
#include <fest/base/time_iterator.tpl.hh>

//#include <fest/numerics/vector_tools.hh>


#include <wave/Problem/Problem.hh>
#include <wave/parameters/ParameterSet.hh>
#include <wave/Grid/TriaGenerator.tpl.hh>

// iDEAL.II includes
#include <ideal.II/base/quadrature_lib.hh>

#include <ideal.II/grid/fixed_tria.hh>

#include <ideal.II/lac/spacetime_vector.hh>


// DEAL.II includes
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/symmetric_tensor.h>


// C++ includes
#include <iostream>
#include <fstream>
#include <sys/stat.h>  // for mkdir



template <int dim>
class ElasticWave_cG : public wave::Problem
{
private:

    dealii::SymmetricTensor<2,dim> get_strain(const dealii::Tensor<2,dim> &grad);

    void setup_dofs_on_slab();
    void assemble_system_on_slab();
    void apply_initial_values(unsigned int start, unsigned int end);
    void apply_final_values();
    void apply_temporal_conditions();
    void solve_on_slab();
    void output_results(const unsigned int refinement_cycle);
    void plot_solution(const unsigned int refinement_cycle);
    void time_steps(const unsigned int refinement_cycle);

    void process_solution(unsigned int cycle);

    void print_convergence_table();

    std::shared_ptr< wave::ParameterSet > parameter_set;

    // struct holding functions
    struct {
        std::shared_ptr< dealii::Function<dim> > inital_value_function;
        std::shared_ptr< dealii::Function<dim> > final_value_function;
        std::shared_ptr< dealii::Function<dim> > boundary_function;
        std::shared_ptr< dealii::TensorFunction<1, dim, double> > rhs_function;
        std::shared_ptr< dealii::Function<dim> > exact_solution;
    } function;

    // struct holding iterators for the slab
    struct {
        idealii::slab::TriaIterator<dim>       tria;
        fest::cg::slab::DoFHandlerIterator<dim> dof;
        idealii::slab::VectorIterator<double>  solution;
    } slab_it;
    
    idealii::spacetime::fixed::Triangulation<dim> triangulation;


    virtual void init_grid();
    virtual void init_functions();

    fest::cg::DoFHandler<dim> dof_handler;

    fest::cg::FiniteElement<dim> fe;
    

    // Space-Time
    dealii::SparseMatrix<double> slab_system_matrix;
    dealii::Vector<double> slab_system_rhs;
    dealii::Vector<double> slab_initial_value;

    idealii::spacetime::Vector<double> solution;

    const dealii::FEValuesExtractors::Vector displacement = 0;
    const dealii::FEValuesExtractors::Vector velocity = dim;

    dealii::types::global_dof_index n_space_u;
    dealii::types::global_dof_index n_space_v;

    dealii::SparsityPattern      sparsity_pattern;
    std::shared_ptr<dealii::AffineConstraints<double>> slab_constraints;

    unsigned int slab;

    dealii::ConvergenceTable convergence_table;
    std::vector<double> L2_error_vals;
    double L2_error;

    fest::RefinementTable refinement_table;


public:
    ElasticWave_cG(unsigned int s = 1, 
                   unsigned int r = 1);
    virtual ~ElasticWave_cG() = default;
    void set_input_parameters(std::shared_ptr< dealii::ParameterHandler > parameter_handler) override;
    void run() override;
};

#endif