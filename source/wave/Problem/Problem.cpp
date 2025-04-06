// PROJECT includes
#include <wave/Problem/Problem.hh>


// DEAL.II includes
#include <deal.II/base/exceptions.h>


// C++ includes


namespace wave {

void Problem::set_input_parameters([[maybe_unused]] std::shared_ptr< dealii::ParameterHandler > parameter_handler) {
    AssertThrow(false, dealii::ExcNotImplemented());
}

void Problem::run() {
    AssertThrow(false, dealii::ExcNotImplemented());
}





}