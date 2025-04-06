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


#ifndef __SpaceTime_dof_handler_tpl_hh
#define __SpaceTime_dof_handler_tpl_hh

// PROJECT includes
#include <fest/grid/SpaceTimeTria.tpl.hh>
#include <fest/fe/fe.tpl.hh>
#include <fest/dofs/Slab_dof_handler.tpl.hh>

// iDEAL.II includes
#include <ideal.II/grid/spacetime_tria.hh>


// DEAL.II includes
#include <deal.II/dofs/dof_handler.h>


// C++ includes
#include <list>
#include <memory>

namespace fest {
namespace cgcg {

template <int dim>
class DoFHandler {

public:
    DoFHandler(fest::cgcg::Triangulation<dim> &tria);

    DoFHandler(const DoFHandler<dim> &other);
    
    void distribute_dofs(fest::cgcg::FiniteElement<dim> fe);

    std::shared_ptr< dealii::DoFHandler<dim> > spatial();
    std::shared_ptr< dealii::DoFHandler<1> > temporal();

    unsigned int n_dofs_spacetime();
    unsigned int n_dofs_space();
    unsigned int n_dofs_time();

    unsigned int dofs_per_cell_time();
    unsigned int dofs_per_cell_space();

    const dealii::IndexSet & locally_owned_dofs();
    
    
    

private:
    std::shared_ptr< dealii::DoFHandler<dim> >  _spatial_dof_handler;
    std::shared_ptr< dealii::DoFHandler<1> >    _temporal_dof_handler;
    dealii::IndexSet _locally_owned_dofs;



};
} // namespace cgcg

namespace cg {

template <int dim>
class DoFHandler
{

public:
     
    DoFHandler(idealii::spacetime::Triangulation<dim> *tria);
  
  
    void
    generate();

    void
    reinit();  
     
    unsigned int
    M();
  
    
    slab::DoFHandlerIterator<dim>
    begin();
      
    slab::DoFHandlerIterator<dim>
    end();
  
protected:
    idealii::spacetime::Triangulation<dim> *_tria;
  
    std::list<slab::DoFHandler<dim>> _dof_handlers;
};

} // namespace cg
} // namespace fest



#endif