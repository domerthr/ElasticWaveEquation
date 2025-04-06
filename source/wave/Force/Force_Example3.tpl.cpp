// PROJECT includes
#include <wave/Force/Force_Example3.tpl.hh>

// DEAL.II includes

// iDEAL.II includes

// C++ includes

namespace wave {
    

template <int dim>
dealii::Tensor<1, dim> ForceExample3<dim>::value([[maybe_unused]] const dealii::Point<dim> &p) const {
    
    dealii::Tensor<1, dim> value;
    const double t = this->get_time();
    value[1] = std::sin(t);
    return value; 
 
}


} // namespace wave 

#include "Force_Example3.inst.in"