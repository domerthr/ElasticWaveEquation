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
#include <fest/base/time_iterator.tpl.hh>


// iDEAL.II includes



// DEAL.II includes



// C++ includes



namespace fest::cg
{
template <int dim>
TimeIteratorCollection<dim>::TimeIteratorCollection()
{
    tria.it_collection  = std::vector<idealii::slab::TriaIterator<dim> *>();
    tria.obj_collection = std::vector<idealii::spacetime::Triangulation<dim> *>();

    dof.it_collection  = std::vector<slab::DoFHandlerIterator<dim> *>();
    dof.obj_collection = std::vector<fest::cg::DoFHandler<dim> *>();

    vector_double.it_collection = std::vector<idealii::slab::VectorIterator<double> *>();

    vector_double.obj_collection = std::vector<idealii::spacetime::Vector<double> *>();

}

template <int dim>
void
TimeIteratorCollection<dim>::add_iterator(idealii::slab::TriaIterator<dim>       *it,
                                          idealii::spacetime::Triangulation<dim> *obj)
{
    tria.it_collection.push_back(it);
    tria.obj_collection.push_back(obj);
}



template <int dim>
void
TimeIteratorCollection<dim>::add_iterator(slab::DoFHandlerIterator<dim> *it,
                                        fest::cg::DoFHandler<dim>    *obj)
{
    dof.it_collection.push_back(it);
    dof.obj_collection.push_back(obj);
}

template <int dim>
void
TimeIteratorCollection<dim>::add_iterator(idealii::slab::VectorIterator<double> *it,
                                          idealii::spacetime::Vector<double>    *obj)
{
    vector_double.it_collection.push_back(it);
    vector_double.obj_collection.push_back(obj);
}

template <int dim>
void
TimeIteratorCollection<dim>::increment()
{
    for (unsigned int i = 0; i < tria.it_collection.size(); i++)
    {
        Assert(*tria.it_collection[i] != tria.obj_collection[i]->end(), dealii::ExcIteratorPastEnd());
        ++(*tria.it_collection[i]);
    }


    for (unsigned int i = 0; i < dof.it_collection.size(); i++)
    {
        Assert(*dof.it_collection[i] != dof.obj_collection[i]->end(),
            dealii::ExcIteratorPastEnd());
        ++(*dof.it_collection[i]);
    }
    for (unsigned int i = 0; i < vector_double.it_collection.size(); i++)
    {
        Assert(*vector_double.it_collection[i] != vector_double.obj_collection[i]->end(), dealii::ExcIteratorPastEnd());
        ++(*vector_double.it_collection[i]);
    }

}

template <int dim>
void
TimeIteratorCollection<dim>::decrement()
{
    for (unsigned int i = 0; i < tria.it_collection.size(); i++)
    {
        Assert(*tria.it_collection[i] !=
                std::prev(tria.obj_collection[i]->begin()),
            dealii::ExcIteratorPastEnd());
        --(*tria.it_collection[i]);
    }


    for (unsigned int i = 0; i < dof.it_collection.size(); i++)
    {
        Assert(*dof.it_collection[i] !=
                std::prev(dof.obj_collection[i]->begin()),
            dealii::ExcIteratorPastEnd());
        --(*dof.it_collection[i]);
    }
    for (unsigned int i = 0; i < vector_double.it_collection.size(); i++)
    {
        Assert(*vector_double.it_collection[i] !=
                std::prev(vector_double.obj_collection[i]->begin()),
            dealii::ExcIteratorPastEnd());
        --(*vector_double.it_collection[i]);
    }

}

template <int dim>
bool
TimeIteratorCollection<dim>::at_end()
{
    bool res = false;
    for (unsigned int i = 0; i < tria.it_collection.size(); i++)
    {
        if (*tria.it_collection[i] == tria.obj_collection[i]->end())
        {
            res = true;
        }
    }
    for (unsigned int i = 0; i < dof.it_collection.size(); i++)
    {
        if (*dof.it_collection[i] == dof.obj_collection[i]->end())
        {
            res = true;
        }
    }
    for (unsigned int i = 0; i < vector_double.it_collection.size(); i++)
    {
        if (*vector_double.it_collection[i] == vector_double.obj_collection[i]->end())
        {
            res = true;
        }
    }
    return res;
}

template <int dim>
bool
TimeIteratorCollection<dim>::before_begin() {
    bool res = false;
    for (unsigned int i = 0; i < tria.it_collection.size(); i++)
        {
        if (*tria.it_collection[i] == std::prev(tria.obj_collection[i]->begin()))
        {
            res = true;
        }
    }

    for (unsigned int i = 0; i < dof.it_collection.size(); i++)
    {
        if (*dof.it_collection[i] == std::prev(dof.obj_collection[i]->begin()))
        {
            res = true;
        }
    }
    for (unsigned int i = 0; i < vector_double.it_collection.size(); i++)
    {
        if (*vector_double.it_collection[i] == std::prev(vector_double.obj_collection[i]->begin()))
        {
            res = true;
        }
    }

    return res;
}
} // namespace idealii

#include "time_iterator.inst.in"