#ifndef __ParameterHandler_hh
#define __ParameterHandler_hh

// PROJECT includes

// DEALI.II includes
#include <deal.II/base/parameter_handler.h>

// C++ includes

namespace wave {

class ParameterHandler : public dealii::ParameterHandler
{
public:
    ParameterHandler();
    virtual ~ParameterHandler() = default;

};

}

#endif