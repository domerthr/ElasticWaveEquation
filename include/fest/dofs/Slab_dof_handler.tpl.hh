// ---------------------------------------------------------------------
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
// ---------------------------------------------------------------------

#ifndef __Slab_dof_handler_tpl_hh
#define __Slab_dof_handler_tpl_hh


// PROJECT includes
#include <fest/fe/fe.tpl.hh>

// iDEAL.II includes
#include <ideal.II/grid/slab_tria.hh>

// DEAL.II includes
#include <deal.II/dofs/dof_handler.h>


// C++ includes
#include <list>
#include <memory>

namespace fest {
namespace cg {
namespace slab 
{

template <int dim>
class DoFHandler
{
public:
 
    DoFHandler(idealii::slab::Triangulation<dim> &tria);


    DoFHandler(const DoFHandler<dim> &other);
    

    std::shared_ptr<dealii::DoFHandler<dim>>
    spatial();

   
    std::shared_ptr<dealii::DoFHandler<1>>
    temporal_trial();

    std::shared_ptr<dealii::DoFHandler<1>>
    temporal_test();

   
    void
    distribute_dofs(fest::cg::FiniteElement<dim> fe);


    unsigned int
    n_dofs_spacetime_trial();

    unsigned int
    n_dofs_spacetime_test();
 

    unsigned int
    n_dofs_space();
    

    unsigned int
    n_dofs_time_trial();

    unsigned int
    n_dofs_time_test();


    unsigned int
    dofs_per_cell_time_trial();

    unsigned int
    dofs_per_cell_time_test();

    
    unsigned int
    dofs_per_cell_space();

    const dealii::IndexSet &
    locally_owned_dofs_trial();

    const dealii::IndexSet &
    locally_owned_dofs_test();

    typename fest::cg::FiniteElement<dim>::support_type
    fe_support_type();

private:
    std::shared_ptr<dealii::DoFHandler<dim>>                _spatial_dof;
    std::shared_ptr<dealii::DoFHandler<1>>                  _temporal_dof_trial;
    std::shared_ptr<dealii::DoFHandler<1>>                  _temporal_dof_test;
    dealii::IndexSet                                        _locally_owned_dofs_trial;
    dealii::IndexSet                                        _locally_owned_dofs_test;
    typename fest::cg::FiniteElement<dim>::support_type     _fe_support_type;
};

 
template <int dim>
using DoFHandlerIterator = typename std::list<DoFHandler<dim>>::iterator;


} // namespace slab
} // namespace cg
} // namespace fest

#endif