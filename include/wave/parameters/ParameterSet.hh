#ifndef __ParameterSet_hh
#define __ParameterSet_hh

// PROJECT includes

// DEALI.II includes
#include <deal.II/base/parameter_handler.h>

// C++ includes
#include <string>

namespace wave {
struct ParameterSet {
    ParameterSet(std::shared_ptr< dealii::ParameterHandler > handler);

    bool junker;
    unsigned int dim;
    unsigned int example;

    unsigned int s;
    unsigned int r;

    // space mesh
    std::string TriaClass;
    std::string TriaOptions;
    unsigned int global_refinement_space;
    bool refine_space;

    // time mesh
    double t0;
    double T;
    unsigned int M;
    unsigned int global_refinement_time;
    bool refine_time;
    unsigned int max_refinements_time;

    // parameter specification
    double lame_coefficient_mu;
    double poisson_ratio_nu;
    double rho;
    double omega_dot;

    


};

}

#endif