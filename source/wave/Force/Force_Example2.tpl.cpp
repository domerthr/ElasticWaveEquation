// PROJECT includes
#include <wave/Force/Force_Example2.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

template <int dim>
dealii::Tensor<1, dim> ForceExample2<dim>::value([[maybe_unused]] const dealii::Point<dim> &p) const {
    
    dealii::Tensor<1, dim> value;
    return value; 
 
}


} // namespace wave 

#include "Force_Example2.inst.in"