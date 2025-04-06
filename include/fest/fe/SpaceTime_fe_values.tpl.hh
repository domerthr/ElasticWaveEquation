#ifndef __SpaceTime_fe_values_tpl_hh
#define __SpaceTime_fe_values_tpl_hh

// PROJECT includes
#include <fest/fe/fe.tpl.hh>
#include <fest/base/quadrature_lib.tpl.hh>

// iDEAL.II includes
#include <ideal.II/base/quadrature_lib.hh>


// DEAL.II includes
#include <deal.II/fe/fe_values.h>


// C++ includes

namespace fest {
namespace cgcg {

template <int dim>
class FEValues 
{

public:
    FEValues(fest::cgcg::FiniteElement<dim> &fe,
             fest::Quadrature<dim> &quad,
             const dealii::UpdateFlags uflags);

    void reinit_space(const typename dealii::TriaIterator< dealii::DoFCellAccessor<dim, dim, false> > &cell_space);
    void reinit_time(const typename dealii::TriaIterator< dealii::DoFCellAccessor<1, 1, false> > &cell_time);

    // get value / gradient
    double shape_value(unsigned int shape_function, unsigned int q_point);

    double shape_dt(unsigned int shape_function, unsigned int q_point);

    typename dealii::FEValuesViews::Scalar<dim>::value_type
    scalar_value(const typename dealii::FEValuesExtractors::Scalar &extractor,
                 unsigned int shape_function,
                 unsigned int q_point);

    typename dealii::FEValuesViews::Scalar<dim>::value_type
    scalar_dt(const typename dealii::FEValuesExtractors::Scalar &extractor,
              unsigned int shape_function,
              unsigned int q_point);

    typename dealii::FEValuesViews::Scalar<dim>::gradient_type
    scalar_space_grad(const typename dealii::FEValuesExtractors::Scalar &extractor,
                      unsigned int shape_function,
                      unsigned int q_point);

    typename dealii::FEValuesViews::Vector<dim>::value_type
    vector_value(const typename dealii::FEValuesExtractors::Vector &extractor,
                 unsigned int shape_function,
                 unsigned int q_point);
    
    typename dealii::FEValuesViews::Vector<dim>::value_type
    vector_dt(const typename dealii::FEValuesExtractors::Vector &extractor,
              unsigned int shape_function,
              unsigned int q_point);
    
    typename dealii::FEValuesViews::Vector<dim>::divergence_type
    vector_divergence(const typename dealii::FEValuesExtractors::Vector &extractor,
                      unsigned int shape_function,
                      unsigned int q_point);


    typename dealii::FEValuesViews::Vector<dim>::gradient_type
    vector_space_grad(const typename dealii::FEValuesExtractors::Vector &extractor,
                      unsigned int shape_function,
                      unsigned int q_point);

    
    double time_quadratur_point(unsigned int q_point);

    dealii::Point<dim> space_quadratur_point(unsigned int q_point);

    double JxW(const unsigned int q_point);

    void get_local_dof_indices(std::vector<dealii::types::global_dof_index> &indices);

    std::shared_ptr< dealii::FEValues<dim> > spatial();
    std::shared_ptr< dealii::FEValues<1> > temporal();

private:
    fest::cgcg::FiniteElement<dim> &_fe;
    fest::Quadrature<dim> &_quad;

    std::shared_ptr< dealii::FEValues<dim> > _fev_space;
    std::shared_ptr< dealii::FEValues<1> > _fev_time;

    std::vector< dealii::types::global_dof_index > local_space_dof_index;
    std::vector< dealii::types::global_dof_index > local_time_dof_index;

    unsigned int       n_dofs_space;
    unsigned int       time_cell_index;
    const unsigned int n_dofs_space_cell;
    const unsigned int n_quads_space;

public:
    unsigned int n_quadature_points;

};

} // namespace cgcg
namespace cg {

template <int dim>
class FEValues
{
public:
    
    FEValues(FiniteElement<dim>    &fe,
             idealii::spacetime::Quadrature<dim>          &quad,
             const dealii::UpdateFlags uflags);

   
    void
    reinit_space(const typename dealii::TriaIterator<
                 dealii::DoFCellAccessor<dim, dim, false>> &cell_space);
    
    void
    reinit_time(const typename dealii::TriaIterator<dealii::DoFCellAccessor<1, 1, false>> &cell_time,
                const typename dealii::TriaIterator<dealii::DoFCellAccessor<1, 1, false>> &cell_time_test);
    
    typename dealii::FEValuesViews::Scalar<dim>::value_type
    scalar_value(const typename dealii::FEValuesExtractors::Scalar &extractor,
                 unsigned int                                      function_no,
                 unsigned int                                      point_no);
    

   
    typename dealii::FEValuesViews::Scalar<dim>::value_type
    scalar_dt(const typename dealii::FEValuesExtractors::Scalar &extractor,
            unsigned int                                       function_no,
            unsigned int                                       point_no);
    
    typename dealii::FEValuesViews::Scalar<dim>::gradient_type
    scalar_space_grad(
    const typename dealii::FEValuesExtractors::Scalar &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no);


   
    typename dealii::FEValuesViews::Vector<dim>::value_type
    vector_value(const typename dealii::FEValuesExtractors::Vector &extractor,
                unsigned int                                       function_no,
                unsigned int                                       point_no);
    
   
    typename dealii::FEValuesViews::Vector<dim>::value_type
    vector_dt(const typename dealii::FEValuesExtractors::Vector &extractor,
              unsigned int                                       function_no,
              unsigned int                                       point_no);

 
    typename dealii::FEValuesViews::Vector<dim>::divergence_type
    vector_divergence(const typename dealii::FEValuesExtractors::Vector &extractor,
                      unsigned int                                       function_no,
                      unsigned int                                       point_no);

   
    typename dealii::FEValuesViews::Vector<dim>::gradient_type
    vector_space_grad(const typename dealii::FEValuesExtractors::Vector &extractor,
                      unsigned int                                       function_no,
                      unsigned int                                       point_no);
    
    typename dealii::FEValuesViews::Vector<dim>::value_type
    vector_value_test(const typename dealii::FEValuesExtractors::Vector &extractor,
                      unsigned int                                       function_no,
                      unsigned int                                       point_no);
    
    
    typename dealii::FEValuesViews::Vector<dim>::value_type
    vector_dt_test(const typename dealii::FEValuesExtractors::Vector &extractor,
                   unsigned int                                       function_no,
                   unsigned int                                       point_no);


    typename dealii::FEValuesViews::Vector<dim>::divergence_type
    vector_divergence_test(const typename dealii::FEValuesExtractors::Vector &extractor,
                           unsigned int                                       function_no,
                           unsigned int                                       point_no);

    
    typename dealii::FEValuesViews::Vector<dim>::gradient_type
    vector_space_grad_test(const typename dealii::FEValuesExtractors::Vector &extractor,
                           unsigned int                                       function_no,
                           unsigned int                                       point_no);

    
 
    double
    time_quadrature_point(unsigned int quadrature_point);


    dealii::Point<dim>
    space_quadrature_point(unsigned int quadrature_point);

   
    double
    JxW(const unsigned int quadrature_point);

    void
    get_local_dof_indices(
    std::vector<dealii::types::global_dof_index> &indices);

 
    std::shared_ptr<dealii::FEValues<dim>>
    spatial();
 
    std::shared_ptr<dealii::FEValues<1>>
    temporal_trial();

    std::shared_ptr<dealii::FEValues<1>>
    temporal_test();

private:
    FiniteElement<dim> &_fe;
    idealii::spacetime::Quadrature<dim>    &_quad;

    std::shared_ptr<dealii::FEValues<dim>> _fev_space;
    std::shared_ptr<dealii::FEValues<1>>   _fev_time_trial;
    std::shared_ptr<dealii::FEValues<1>>   _fev_time_test;

    std::vector<dealii::types::global_dof_index> local_space_dof_index;
    std::vector<dealii::types::global_dof_index> local_time_dof_index_trial;
    std::vector<dealii::types::global_dof_index> local_time_dof_index_test;

    unsigned int       n_dofs_space;
    unsigned int       time_cell_index_trial;
    unsigned int       time_cell_index_test;
    const unsigned int n_dofs_space_cell;
    const unsigned int n_quads_space;

public:
    
    unsigned int n_quadrature_points;
};



} // namespace cg

} // namespace fest

#endif