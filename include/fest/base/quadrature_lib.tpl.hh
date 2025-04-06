#ifndef __quadrature_lib_tpl_hh
#define __quadrature_lib_tpl_hh

// PROJECT includes
#include <fest/base/SpaceTime_quadrature.tpl.hh>


// DEAL.II includes
#include <deal.II/base/quadrature_lib.h>


// C++ includes

namespace fest {

template <int dim>
class QGauss : public fest::Quadrature<dim>
{

public:
    QGauss(unsigned int n_space, unsigned int n_time);
    

};




}

#endif