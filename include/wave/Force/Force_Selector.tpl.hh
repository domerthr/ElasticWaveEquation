#ifndef __Force_Selector_tpl_hh
#define __Force_Selector_tpl_hh

// DEALI.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

// PROJECT includes


// C++ includes
#include <string>

namespace wave {
namespace force {

template <int dim>
class Selector {

public:
    Selector() = default;
    virtual ~Selector() = default;

    virtual void create_function(const unsigned int example,
                                 std::shared_ptr < dealii::TensorFunction<1, dim, double>> &function) const;
};

}}

#endif