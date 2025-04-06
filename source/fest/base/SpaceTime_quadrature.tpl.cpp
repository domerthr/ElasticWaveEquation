// PROJECT includes
#include <fest/base/SpaceTime_quadrature.tpl.hh>

// DEAL.II includes



// C++ includes

namespace fest {

template <int dim>
Quadrature<dim>::Quadrature(std::shared_ptr< dealii::Quadrature<dim> > quad_space,
                            std::shared_ptr< dealii::Quadrature<1> > quad_time) {
    _quad_space = quad_space;
    _quad_time = quad_time;                            

}

template <int dim>    
std::shared_ptr< dealii::Quadrature<dim> > Quadrature<dim>::spatial() {
    return _quad_space;
}

template <int dim>
std::shared_ptr< dealii::Quadrature<1> > Quadrature<dim>::temporal() {
    return _quad_time;
}

}

#include "SpaceTime_quadrature.inst.in"
