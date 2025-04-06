// PROJECT includes
#include <fest/grid/SpaceTimeTria.tpl.hh>


// DEAL.II includes



// C++ includes

namespace fest {
namespace cgcg {

template <int dim>
Triangulation<dim>::Triangulation() {
    _temporal_tria = std::make_shared< dealii::Triangulation<1> >();
    _spatial_tria = std::make_shared< dealii::Triangulation<dim> >();
}

template <int dim>
void Triangulation<dim>::generate(double startpoint, double endpoint) 
{
    _startpoint = startpoint;
    _endpoint = endpoint;
    dealii::GridGenerator::hyper_cube(*_temporal_tria, _startpoint, _endpoint);

}

template <int dim>
Triangulation<dim>::Triangulation(const Triangulation &other)
    :   _startpoint(other._startpoint)
    ,   _endpoint(other._endpoint)
{
    Assert(other._spatial_tria.use_count(), dealii::ExcNotInitialized());
    _spatial_tria = other._spatial_tria;
    Assert(other._temporal_tria.use_count(), dealii::ExcNotInitialized());
    _temporal_tria = other._temporal_tria;
}

template <int dim>
std::shared_ptr< dealii::Triangulation<dim> > Triangulation<dim>::spatial() {
    Assert(_spatial_tria.use_count(), dealii::ExcNotInitialized());
    return _spatial_tria;
}

template <int dim>
std::shared_ptr< dealii::Triangulation<1> > Triangulation<dim>::temporal() {
    Assert(_temporal_tria.use_count(), dealii::ExcNotInitialized());
    return _temporal_tria;
}

template <int dim>
double Triangulation<dim>::startpoint() {
    return _startpoint;
}

template <int dim>
double Triangulation<dim>::endpoint() {
    return _endpoint;
}


template <int dim>
void Triangulation<dim>::refine_global(const unsigned int n, const unsigned int m) {
    
    _spatial_tria->refine_global(n);

    _temporal_tria->refine_global(m);

    
}

template <int dim>
unsigned int Triangulation<dim>::n_active_time_cells() {
    return _temporal_tria->n_active_cells();
}

template <int dim>
unsigned int Triangulation<dim>::n_active_space_cells() {
    return _spatial_tria->n_active_cells();
}

template <int dim>
unsigned int Triangulation<dim>::n_active_cells() {
    return _temporal_tria->n_active_cells() * _spatial_tria->n_active_cells();
}



} // namespace cgcg
} // namespace fest

#include "SpaceTimeTria.inst.in"