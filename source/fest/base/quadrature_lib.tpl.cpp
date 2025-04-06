// PROJECT includes
#include <fest/base/quadrature_lib.tpl.hh>


// DEAL.II includes
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>


// C++ includes
#include <cmath>
#include <vector>

namespace fest {

template <int dim>
QGauss<dim>::QGauss(unsigned int n_space, unsigned int n_time) 
    :   fest::Quadrature<dim>(std::make_shared< dealii::QGauss<dim> >(n_space),
                              std::make_shared< dealii::QGauss<1> >(n_time))
{ }


}

#include "quadrature_lib.inst.in"