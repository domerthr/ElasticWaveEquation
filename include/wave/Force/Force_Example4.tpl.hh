#ifndef __Force_Example4_tpl_hh
#define __Force_Example4_tpl_hh

// PROJECT includes

// DEAL.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/point.h>


// C++ includes

namespace wave {
    
template <int dim>
class ForceExample4 : public dealii::TensorFunction<1, dim>
{

public:
    ForceExample4() 
            : dealii::TensorFunction<1, dim> () {};

    virtual ~ForceExample4() = default;
    
    virtual dealii::Tensor<1, dim> value(const dealii::Point<dim> &p) const override;
        

};
} // namespace wave

#endif