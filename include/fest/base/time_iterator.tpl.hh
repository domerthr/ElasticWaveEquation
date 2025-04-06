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

#ifndef __time_iterator_tpl_hh
#define __time_iterator_tpl_hh

// PROJECT includes
#include <fest/dofs/SpaceTime_dof_handler.tpl.hh>


// iDEAL.II includes
#include <ideal.II/grid/spacetime_tria.hh>
#include <ideal.II/lac/spacetime_vector.hh>


// DEAL.II includes



// C++ includes
#include <list>
#include <memory>


namespace fest::cg
{


template <int dim>
class TimeIteratorCollection
{
public:

    TimeIteratorCollection();

    void
    add_iterator(idealii::slab::TriaIterator<dim>        *it,
                 idealii::spacetime::Triangulation<dim>  *collection);


    void
    add_iterator(slab::DoFHandlerIterator<dim> *it,
                 fest::cg::DoFHandler<dim>     *collection);


    void
    add_iterator(idealii::slab::VectorIterator<double>  *it,
                 idealii::spacetime::Vector<double>     *collection);


    void
    increment();
 
    void
    decrement();

    bool
    at_end();

    bool
    before_begin();

    private:
    struct
    {
        std::vector<idealii::slab::TriaIterator<dim> *>        it_collection;
        std::vector<idealii::spacetime::Triangulation<dim> *>  obj_collection;
    } tria;


    struct
    {
        std::vector<slab::DoFHandlerIterator<dim> *> it_collection;
        std::vector<fest::cg::DoFHandler<dim> *>     obj_collection;
    } dof;
    struct
    {
        std::vector<idealii::slab::VectorIterator<double> *>  it_collection;
        std::vector<idealii::spacetime::Vector<double> *>     obj_collection;
    } vector_double;
};

} // namespace fest::cg

#endif 