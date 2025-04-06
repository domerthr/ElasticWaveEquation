#ifndef __Force_Example3_tpl_hh
#define __Force_Example3_tpl_hh

// PROJECT includes

// DEAL.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/point.h>


// C++ includes

namespace wave {
    
template <int dim>
class ForceExample3 : public dealii::TensorFunction<1, dim>
{


public:
    ForceExample3() 
            : dealii::TensorFunction<1, dim> () {};

    virtual ~ForceExample3() = default;
    
    virtual dealii::Tensor<1, dim> value(const dealii::Point<dim> &p) const override;
        

};
} // namespace wave

#endif