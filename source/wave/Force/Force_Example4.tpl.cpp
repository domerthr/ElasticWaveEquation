// PROJECT includes
#include <wave/Force/Force_Example4.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

template <int dim>
dealii::Tensor<1, dim> ForceExample4<dim>::value(const dealii::Point<dim> &p) const {
    
    dealii::Tensor<1, dim> value;
    const double t = this->get_time();
    value[0] = (M_PI * M_PI - 1) * std::sin(t) * std::sin(M_PI * p[0]);
    return value; 
 
}


} // namespace wave 

#include "Force_Example4.inst.in"