#ifndef __wave_Problem_hh
#define __wave_Problem_hh

// PROJECT includes



// DEAL.II includes
#include <deal.II/base/parameter_handler.h>


// C++ includes
#include <memory>

namespace wave {

class Problem {
public:
    Problem() = default;
    virtual ~Problem() = default;

    virtual void set_input_parameters(std::shared_ptr< dealii::ParameterHandler > parameter_handler);
    virtual void run();

};

} // namespace


#endif