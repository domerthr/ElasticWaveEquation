// PROJECT includes
#include <wave/Force/Force_Example1.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

template <int dim>
dealii::Tensor<1, dim> ForceExample1<dim>::value(const dealii::Point<dim> &p) const {
    
    dealii::Tensor<1, dim> value;
    const double t = this->get_time();
    value[0] = (2. - p[0] * (1. - p[0])) * std::sin(t);
    return value; 
 
}


} // namespace wave 

#include "Force_Example1.inst.in"