// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the ideal.II authors
//
// This file is based on the ideal.II library.
//
// The ideal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of ideal.II.
//
// ---------------------------------------------------------------------



// PROJECT includes
#include <fest/dofs/Slab_dof_handler.tpl.hh>


// DEAL.II includes



// C++ includes


namespace fest {
namespace cg {
namespace slab 
{

template <int dim>
DoFHandler<dim>::DoFHandler(idealii::slab::Triangulation<dim> &tria) {
    Assert(tria.spatial().use_count(), dealii::ExcNotInitialized());
    Assert(tria.temporal().use_count(), dealii::ExcNotInitialized());
    _spatial_dof  = std::make_shared<dealii::DoFHandler<dim>>(*tria.spatial());
    _temporal_dof_trial = std::make_shared<dealii::DoFHandler<1>>(*tria.temporal());
    _temporal_dof_test = std::make_shared<dealii::DoFHandler<1>>(*tria.temporal());
    _locally_owned_dofs_trial = dealii::IndexSet();
    _locally_owned_dofs_test = dealii::IndexSet();
}

template <int dim>
DoFHandler<dim>::DoFHandler(const DoFHandler<dim> &other) {
    Assert(other._spatial_dof.use_count(), dealii::ExcNotInitialized());
    _spatial_dof = other._spatial_dof;
    Assert(other._temporal_dof_trial.use_count(), dealii::ExcNotInitialized());
    Assert(other._temporal_dof_test.use_count(), dealii::ExcNotInitialized());
    _temporal_dof_trial       = other._temporal_dof_trial;
    _temporal_dof_test        = other._temporal_dof_test;
    _locally_owned_dofs_trial = other._locally_owned_dofs_trial;
    _locally_owned_dofs_test = other._locally_owned_dofs_test;
    _fe_support_type    = other._fe_support_type;
}

template <int dim>
std::shared_ptr<dealii::DoFHandler<dim>> DoFHandler<dim>::spatial() {
    return _spatial_dof;
}

template <int dim> 
std::shared_ptr<dealii::DoFHandler<1>> DoFHandler<dim>::temporal_trial() {
    return _temporal_dof_trial;
}

template <int dim> 
std::shared_ptr<dealii::DoFHandler<1>> DoFHandler<dim>::temporal_test() {
    return _temporal_dof_test;
}

template <int dim>
void DoFHandler<dim>::distribute_dofs(fest::cg::FiniteElement<dim> fe) {

    _fe_support_type = fe.type();
    _spatial_dof->distribute_dofs(*fe.spatial());
    _temporal_dof_trial->distribute_dofs(*fe.temporal_trial());

    _temporal_dof_test->distribute_dofs(*fe.temporal_test());
    dealii::IndexSet space_owned_dofs = _spatial_dof->locally_owned_dofs();

    _locally_owned_dofs_trial.clear();
    _locally_owned_dofs_trial.set_size(space_owned_dofs.size() * _temporal_dof_trial->n_dofs());

    for (dealii::types::global_dof_index time_dof_index = 0;
         time_dof_index < _temporal_dof_trial->n_dofs();
         time_dof_index++) {

        _locally_owned_dofs_trial.add_indices(space_owned_dofs,
                                              time_dof_index *
                                              _spatial_dof->n_dofs() // offset
        );
    }

    _locally_owned_dofs_test.clear();
    _locally_owned_dofs_test.set_size(space_owned_dofs.size() * _temporal_dof_test->n_dofs());

    for (dealii::types::global_dof_index time_dof_index = 0;
        time_dof_index < _temporal_dof_test->n_dofs();
        time_dof_index++) {

       _locally_owned_dofs_test.add_indices(space_owned_dofs,
                                            time_dof_index *
                                            _spatial_dof->n_dofs() // offset
       );
    }
    
}

template <int dim>
unsigned int DoFHandler<dim>::n_dofs_spacetime_trial() {
    return _spatial_dof->n_dofs() * _temporal_dof_trial->n_dofs();
}

template <int dim>
unsigned int DoFHandler<dim>::n_dofs_spacetime_test() {
    return _spatial_dof->n_dofs() * _temporal_dof_test->n_dofs();
}

template <int dim>
unsigned int DoFHandler<dim>::n_dofs_space() {

    return _spatial_dof->n_dofs();
}


template <int dim>
unsigned int DoFHandler<dim>::n_dofs_time_trial() {
    
    return _temporal_dof_trial->n_dofs();
}

template <int dim>
unsigned int DoFHandler<dim>::n_dofs_time_test() {
    
    return _temporal_dof_test->n_dofs();
}

template <int dim>
unsigned int DoFHandler<dim>::dofs_per_cell_time_trial() {
    
    return _temporal_dof_trial->get_fe().dofs_per_cell;
}

template <int dim>
unsigned int DoFHandler<dim>::dofs_per_cell_time_test() {
    
    return _temporal_dof_test->get_fe().dofs_per_cell;
}

template <int dim>
unsigned int DoFHandler<dim>::dofs_per_cell_space() {
    
    return _spatial_dof->get_fe().dofs_per_cell;
}

template <int dim>
const dealii::IndexSet &
DoFHandler<dim>::locally_owned_dofs_trial() {
    return _locally_owned_dofs_trial;
}

template <int dim>
const dealii::IndexSet &
DoFHandler<dim>::locally_owned_dofs_test() {
    return _locally_owned_dofs_test;
}

template <int dim>
typename fest::cg::FiniteElement<dim>::support_type DoFHandler<dim>::fe_support_type()
{
    return _fe_support_type;
}


} // namespace slab
} // namespace cg
} // namespace fest


#include "Slab_dof_handler.inst.in"