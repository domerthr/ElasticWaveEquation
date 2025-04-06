// PROJECT includes
#include <fest/fe/SpaceTime_fe_values.tpl.hh>


// DEAL.II includes



// C++ includes

namespace fest {
namespace cgcg {

template <int dim>
FEValues<dim>::FEValues(fest::cgcg::FiniteElement<dim> &fe,
                        fest::Quadrature<dim> &quad,
                        const dealii::UpdateFlags uflags)
    :   _fe(fe)
    ,   _quad(quad)
    ,   _fev_space(std::make_shared< dealii::FEValues<dim> > (*fe.spatial(), *quad.spatial(), uflags))
    ,   _fev_time(std::make_shared< dealii::FEValues<1> > (*fe.temporal(), *quad.temporal(), uflags))
    ,   local_space_dof_index(fe.spatial()->dofs_per_cell)
    ,   local_time_dof_index(fe.temporal()->dofs_per_cell)
    ,   n_dofs_space(0)
    ,   time_cell_index(0)
    ,   n_dofs_space_cell(_fe.spatial()->dofs_per_cell)
    ,   n_quads_space(_fev_space->n_quadrature_points)
    ,   n_quadature_points(_fev_space->n_quadrature_points * _fev_time->n_quadrature_points)                 
{ }

template <int dim>
void FEValues<dim>::reinit_space(const typename dealii::TriaIterator< dealii::DoFCellAccessor<dim, dim, false> > &cell_space) {
    _fev_space->reinit(cell_space);
    cell_space->get_dof_indices(local_space_dof_index);
    n_dofs_space = cell_space->get_dof_handler().n_dofs();
}

template <int dim>
void FEValues<dim>::reinit_time(const typename dealii::TriaIterator< dealii::DoFCellAccessor<1, 1, false> > &cell_time) {
    _fev_time->reinit(cell_time);
    cell_time->get_dof_indices(local_time_dof_index);
    time_cell_index = cell_time->index();
}

template <int dim>
double FEValues<dim>::shape_value(unsigned int shape_function, unsigned int q_point) {
    return  _fev_space->shape_value(shape_function % n_dofs_space_cell, q_point % n_quads_space) * 
            _fev_time->shape_value(shape_function / n_dofs_space_cell, q_point / n_quads_space);
}

template <int dim>
double FEValues<dim>::shape_dt(unsigned int shape_function, unsigned int q_point) {
    return _fev_space->shape_value(shape_function % n_dofs_space_cell, q_point % n_quads_space) *
           _fev_time->shape_grad(shape_function / n_dofs_space_cell, q_point / n_quads_space)[0];
}

template <int dim>
typename dealii::FEValuesViews::Scalar<dim>::value_type
FEValues<dim>::scalar_value(const typename dealii::FEValuesExtractors::Scalar &extractor,
                            unsigned int shape_function,
                            unsigned int q_point) {
    
    return (*_fev_space)[extractor].value(shape_function % n_dofs_space_cell, q_point % n_quads_space) * 
           _fev_time->shape_value(shape_function / n_dofs_space_cell, q_point / n_quads_space);
}

template <int dim>
typename dealii::FEValuesViews::Scalar<dim>::value_type
FEValues<dim>::scalar_dt(const typename dealii::FEValuesExtractors::Scalar &extractor,
                         unsigned int shape_function,
                         unsigned int q_point) {
    
    return (*_fev_space)[extractor].value(shape_function % n_dofs_space_cell, q_point % n_quads_space) * 
           _fev_time->shape_grad(shape_function / n_dofs_space_cell, q_point / n_quads_space)[0];
}

template <int dim>
typename dealii::FEValuesViews::Scalar<dim>::gradient_type
FEValues<dim>::scalar_space_grad(const typename dealii::FEValuesExtractors::Scalar &extractor,
                                 unsigned int shape_function,
                                 unsigned int q_point) {
    
    return (*_fev_space)[extractor].gradient(shape_function % n_dofs_space_cell, q_point % n_quads_space) * 
           _fev_time->shape_value(shape_function / n_dofs_space_cell, q_point / n_quads_space);
}


template <int dim>
typename dealii::FEValuesViews::Vector<dim>::value_type
FEValues<dim>::vector_value(const typename dealii::FEValuesExtractors::Vector &extractor,
                            unsigned int shape_function,
                            unsigned int q_point)
{
    return (*_fev_space)[extractor].value(shape_function % n_dofs_space_cell, q_point % n_quads_space) *
           _fev_time->shape_value(shape_function / n_dofs_space_cell, q_point / n_quads_space);
}

template <int dim>
typename dealii::FEValuesViews::Vector<dim>::value_type
FEValues<dim>::vector_dt(const typename dealii::FEValuesExtractors::Vector &extractor,
                         unsigned int shape_function,
                         unsigned int q_point) {

    return (*_fev_space)[extractor].value(shape_function % n_dofs_space_cell, q_point % n_quads_space) *
            _fev_time->shape_grad(shape_function / n_dofs_space_cell, q_point / n_quads_space)[0];
}

template <int dim>
typename dealii::FEValuesViews::Vector<dim>::divergence_type
FEValues<dim>::vector_divergence(const typename dealii::FEValuesExtractors::Vector &extractor,
                                 unsigned int shape_function,
                                 unsigned int q_point) {
    return (*_fev_space)[extractor].divergence(shape_function % n_dofs_space_cell, q_point % n_quads_space) *
            _fev_time->shape_value(shape_function / n_dofs_space_cell, q_point / n_quads_space);
}

template <int dim>
typename dealii::FEValuesViews::Vector<dim>::gradient_type
FEValues<dim>::vector_space_grad(const typename dealii::FEValuesExtractors::Vector &extractor,
                                 unsigned int shape_function,
                                 unsigned int q_point) {
    return (*_fev_space)[extractor].gradient(shape_function % n_dofs_space_cell, q_point % n_quads_space) *
            _fev_time->shape_value(shape_function / n_dofs_space_cell, q_point / n_quads_space);
}




template <int dim>
double FEValues<dim>::time_quadratur_point(unsigned int q_point) {
    return _fev_time->quadrature_point(q_point / n_quads_space)[0];
}


template <int dim>
dealii::Point<dim> FEValues<dim>::space_quadratur_point(unsigned int q_point) {
    return _fev_space->quadrature_point(q_point % n_quads_space);
}

template <int dim>
double FEValues<dim>::JxW(const unsigned int q_point) {

    return _fev_space->JxW(q_point % n_quads_space) * _fev_time->JxW(q_point / n_quads_space);

}

template <int dim>
void FEValues<dim>::get_local_dof_indices(std::vector<dealii::types::global_dof_index> &indices) {


    for (unsigned int i = 0; i < _fe.dofs_per_cell; i++) {
        indices[i + time_cell_index * _fe.dofs_per_cell] =
          local_space_dof_index[i % n_dofs_space_cell] +
          local_time_dof_index[i / n_dofs_space_cell] * n_dofs_space;
        
    }

}

template <int dim>
std::shared_ptr< dealii::FEValues<dim> > FEValues<dim>::spatial() {
    return _fev_space;
}

template <int dim>
std::shared_ptr< dealii::FEValues<1> > FEValues<dim>::temporal() {
    return _fev_time;
}

} // namespace cgcg

namespace cg {

template <int dim>
FEValues<dim>::FEValues(FiniteElement<dim>    &fe,
                        idealii::spacetime::Quadrature<dim>          &quad,
                        const dealii::UpdateFlags uflags)
    :   _fe(fe)
    ,   _quad(quad)
    ,   _fev_space(std::make_shared<dealii::FEValues<dim>>(*fe.spatial(),
                                                           *quad.spatial(),
                                                            uflags))
    ,   _fev_time_trial(std::make_shared<dealii::FEValues<1>>(*fe.temporal_trial(),
                                                              *quad.temporal(),
                                                              uflags))
    ,   _fev_time_test(std::make_shared<dealii::FEValues<1>>(*fe.temporal_test(),
                                                             *quad.temporal(),
                                                             uflags))
    ,   local_space_dof_index(fe.spatial()->dofs_per_cell)
    ,   local_time_dof_index_trial(fe.temporal_trial()->dofs_per_cell)
    ,   local_time_dof_index_test(fe.temporal_test()->dofs_per_cell)
    ,   n_dofs_space(0)
    ,   time_cell_index_trial(0)
    ,   time_cell_index_test(0)
    ,   n_dofs_space_cell(_fe.spatial()->dofs_per_cell)
    ,   n_quads_space(_fev_space->n_quadrature_points)
    ,   n_quadrature_points(_fev_space->n_quadrature_points * _fev_time_trial->n_quadrature_points)
{}

template <int dim>
void FEValues<dim>::reinit_space(const typename dealii::TriaIterator< dealii::DoFCellAccessor<dim, dim, false> > &cell_space) {
    _fev_space->reinit(cell_space);
    cell_space->get_dof_indices(local_space_dof_index);
    n_dofs_space = cell_space->get_dof_handler().n_dofs();
}

template <int dim>
void
FEValues<dim>::reinit_time(
    const typename dealii::TriaIterator<dealii::DoFCellAccessor<1, 1, false>>
    &cell_time,
    const typename dealii::TriaIterator<dealii::DoFCellAccessor<1, 1, false>>
    &cell_time_test)
{
    _fev_time_trial->reinit(cell_time);
    _fev_time_test->reinit(cell_time_test);
    cell_time->get_dof_indices(local_time_dof_index_trial);
    cell_time_test->get_dof_indices(local_time_dof_index_test);
    time_cell_index_trial = cell_time->index();
    time_cell_index_test = cell_time_test->index();
}


template <int dim>
typename dealii::FEValuesViews::Scalar<dim>::value_type
FEValues<dim>::scalar_value(
    const typename dealii::FEValuesExtractors::Scalar &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].value(function_no % n_dofs_space_cell,
                                        point_no % n_quads_space) *
            _fev_time_trial->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}

template <int dim>
typename dealii::FEValuesViews::Scalar<dim>::value_type
FEValues<dim>::scalar_dt(
    const typename dealii::FEValuesExtractors::Scalar &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].value(function_no % n_dofs_space_cell,
                                        point_no % n_quads_space) *
            _fev_time_trial->shape_grad(function_no / n_dofs_space_cell,
                                point_no / n_quads_space)[0];
}

template <int dim>
typename dealii::FEValuesViews::Scalar<dim>::gradient_type
FEValues<dim>::scalar_space_grad(
    const typename dealii::FEValuesExtractors::Scalar &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].gradient(function_no % n_dofs_space_cell,
                                            point_no % n_quads_space) *
            _fev_time_trial->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}


template <int dim>
typename dealii::FEValuesViews::Vector<dim>::value_type
FEValues<dim>::vector_value(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].value(function_no % n_dofs_space_cell,
                                        point_no % n_quads_space) *
            _fev_time_trial->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}



template <int dim>
typename dealii::FEValuesViews::Vector<dim>::value_type
FEValues<dim>::vector_dt(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].value(function_no % n_dofs_space_cell,
                                        point_no % n_quads_space) *
            _fev_time_trial->shape_grad(function_no / n_dofs_space_cell,
                                point_no / n_quads_space)[0];
}

template <int dim>
typename dealii::FEValuesViews::Vector<dim>::divergence_type
FEValues<dim>::vector_divergence(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].divergence(function_no % n_dofs_space_cell,
                                                point_no % n_quads_space) *
            _fev_time_trial->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}



template <int dim>
typename dealii::FEValuesViews::Vector<dim>::gradient_type
FEValues<dim>::vector_space_grad(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].gradient(function_no % n_dofs_space_cell,
                                            point_no % n_quads_space) *
            _fev_time_trial->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}


template <int dim>
typename dealii::FEValuesViews::Vector<dim>::value_type
FEValues<dim>::vector_value_test(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].value(function_no % n_dofs_space_cell,
                                        point_no % n_quads_space) *
            _fev_time_test->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}



template <int dim>
typename dealii::FEValuesViews::Vector<dim>::value_type
FEValues<dim>::vector_dt_test(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].value(function_no % n_dofs_space_cell,
                                        point_no % n_quads_space) *
            _fev_time_test->shape_grad(function_no / n_dofs_space_cell,
                                point_no / n_quads_space)[0];
}

template <int dim>
typename dealii::FEValuesViews::Vector<dim>::divergence_type
FEValues<dim>::vector_divergence_test(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].divergence(function_no % n_dofs_space_cell,
                                                point_no % n_quads_space) *
            _fev_time_test->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}



template <int dim>
typename dealii::FEValuesViews::Vector<dim>::gradient_type
FEValues<dim>::vector_space_grad_test(
    const typename dealii::FEValuesExtractors::Vector &extractor,
    unsigned int                                       function_no,
    unsigned int                                       point_no)
{
    return (*_fev_space)[extractor].gradient(function_no % n_dofs_space_cell,
                                            point_no % n_quads_space) *
            _fev_time_test->shape_value(function_no / n_dofs_space_cell,
                                point_no / n_quads_space);
}




template <int dim>
double
FEValues<dim>::time_quadrature_point(unsigned int quadrature_point)
{
    return _fev_time_trial->quadrature_point(quadrature_point / n_quads_space)[0];
}

template <int dim>
dealii::Point<dim>
FEValues<dim>::space_quadrature_point(unsigned int quadrature_point)
{
    return _fev_space->quadrature_point(quadrature_point % n_quads_space);
}

template <int dim>
double
FEValues<dim>::JxW(unsigned int quadrature_point)
{
    return _fev_space->JxW(quadrature_point % n_quads_space) *
           _fev_time_trial->JxW(quadrature_point / n_quads_space);
}

template <int dim>
void
FEValues<dim>::get_local_dof_indices(
    std::vector<dealii::types::global_dof_index> &indices)
{
    for (unsigned int i = 0; i < _fe.dofs_per_cell_trial; i++)
    {
        indices[i + time_cell_index_trial * _fe.dofs_per_cell_trial] =
        local_space_dof_index[i % n_dofs_space_cell] +
        local_time_dof_index_trial[i / n_dofs_space_cell] * n_dofs_space;
    }
}

template <int dim>
std::shared_ptr<dealii::FEValues<dim>>
FEValues<dim>::spatial()
{
    return _fev_space;
}

template <int dim>
std::shared_ptr<dealii::FEValues<1>>
FEValues<dim>::temporal_trial()
{
    return _fev_time_trial;
}

template <int dim>
std::shared_ptr<dealii::FEValues<1>>
FEValues<dim>::temporal_test()
{
    return _fev_time_test;
}


} // naemspace cg
} // namespace fest

#include "SpaceTime_fe_values.inst.in"