#ifndef __FinalValue_Selector_tpl_hh
#define __FinalValue_Selector_tpl_hh

// DEALI.II includes
#include <deal.II/base/function.h>

// PROJECT includes


// C++ includes
#include <string>

namespace wave {
namespace final_value {

template <int dim>
class Selector {

public:
    Selector() = default;
    virtual ~Selector() = default;

    virtual void create_function(const unsigned int example,
                                 std::shared_ptr < dealii::Function<dim>> &function,
                                 double T) const;
};

}}

#endif