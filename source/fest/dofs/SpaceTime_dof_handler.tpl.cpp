// PROJECT includes
#include <fest/dofs/SpaceTime_dof_handler.tpl.hh>


// DEAL.II includes



// C++ includes

namespace fest {
namespace cgcg {

template <int dim>
DoFHandler<dim>::DoFHandler(fest::cgcg::Triangulation<dim> &tria) {
    Assert(tria.spatial().use_count(), dealii::ExcNotInitialized());
    Assert(tria.temporal().use_count(), dealii::ExcNotInitialized());
    _spatial_dof_handler = std::make_shared< dealii::DoFHandler<dim> >(*tria.spatial());
    _temporal_dof_handler = std::make_shared< dealii::DoFHandler<1> >(*tria.temporal());
    _locally_owned_dofs = dealii::IndexSet();
}

template <int dim>
DoFHandler<dim>::DoFHandler(const DoFHandler<dim> &other) {

    Assert(other._spatial_dof_handler.use_count(), dealii::ExcNotInitialized());
    _spatial_dof_handler = other._spatial_dof_handler;
    Assert(other._temporal_dof_handler.use_count(), dealii::ExcNotInitialized());
    _temporal_dof_handler       = other._temporal_dof_handler;
    _locally_owned_dofs = other._locally_owned_dofs;

}



template <int dim>
void DoFHandler<dim>::distribute_dofs(fest::cgcg::FiniteElement<dim> fe) {
    _spatial_dof_handler->distribute_dofs(*fe.spatial());
    _temporal_dof_handler->distribute_dofs(*fe.temporal());

    dealii::IndexSet space_owned_dofs = _spatial_dof_handler->locally_owned_dofs();
    _locally_owned_dofs.clear();
    _locally_owned_dofs.set_size(space_owned_dofs.size() *
                                 _temporal_dof_handler->n_dofs());
    

    for (dealii::types::global_dof_index time_dof_index = 0; time_dof_index < _temporal_dof_handler->n_dofs(); time_dof_index++) {
        _locally_owned_dofs.add_indices(space_owned_dofs,
                                        time_dof_index *
                                        _spatial_dof_handler->n_dofs() // offset
        );
    }

}

template <int dim>
std::shared_ptr< dealii::DoFHandler<dim> > DoFHandler<dim>::spatial() {
    return _spatial_dof_handler;
}

template <int dim>
std::shared_ptr< dealii::DoFHandler<1> > DoFHandler<dim>::temporal() {
    return _temporal_dof_handler;
}

template <int dim>
unsigned int DoFHandler<dim>::n_dofs_spacetime() {
    return _spatial_dof_handler->n_dofs() * _temporal_dof_handler->n_dofs();
}

template <int dim>
unsigned int DoFHandler<dim>::n_dofs_space() {
    return _spatial_dof_handler->n_dofs();
}

template <int dim>
unsigned int DoFHandler<dim>::n_dofs_time() {
    return _temporal_dof_handler->n_dofs();
}

template <int dim>
unsigned int DoFHandler<dim>::dofs_per_cell_time() {
    return _temporal_dof_handler->get_fe().dofs_per_cell;
}

template <int dim>
unsigned int DoFHandler<dim>::dofs_per_cell_space() {
    return _spatial_dof_handler->get_fe().dofs_per_cell;
}

template <int dim>
const dealii::IndexSet & DoFHandler<dim>::locally_owned_dofs() {
    return _locally_owned_dofs;
}


} // namespace cgcg


namespace cg {

template <int dim> 
DoFHandler<dim>::DoFHandler(idealii::spacetime::Triangulation<dim> *tria) {
    
    _tria = tria;

    _dof_handlers = std::list<fest::cg::slab::DoFHandler<dim>>();
}

template <int dim>
unsigned int DoFHandler<dim>::M() {
    return _dof_handlers.size();
}

template <int dim>
void DoFHandler<dim>::generate() {

    if (_tria != nullptr) {

        idealii::slab::TriaIterator<dim> tria_it  = this->_tria->begin();
        idealii::slab::TriaIterator<dim> tria_end = this->_tria->end();

        for (; tria_it != tria_end; ++tria_it) {
            this->_dof_handlers.push_back(
            fest::cg::slab::DoFHandler<dim>(*tria_it));
        }
    }

    else {
        Assert(false, dealii::ExcInternalError());
    }
}

template <int dim>
void DoFHandler<dim>::reinit() {
    if (_tria != nullptr) {

        this->_dof_handlers.clear();
        idealii::slab::TriaIterator<dim> tria_it  = this->_tria->begin();
        idealii::slab::TriaIterator<dim> tria_end = this->_tria->end();
        for (; tria_it != tria_end; ++tria_it) {
            this->_dof_handlers.push_back(
                fest::cg::slab::DoFHandler<dim>(*tria_it));
        }
    }
    else {
        Assert(false, dealii::ExcInternalError());
    }
}

template <int dim>
slab::DoFHandlerIterator<dim> DoFHandler<dim>::begin() {

    return _dof_handlers.begin();
}

template <int dim>
slab::DoFHandlerIterator<dim> DoFHandler<dim>::end() {

    return _dof_handlers.end();
}

}


}


#include "SpaceTime_dof_handler.inst.in"