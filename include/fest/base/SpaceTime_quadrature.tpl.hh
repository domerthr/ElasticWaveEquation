#ifndef __SpaceTime_quadrature_tpl_hh
#define __SpaceTime_quadrature_tpl_hh

// PROJECT includes


// DEAL.II includes
#include <deal.II/base/quadrature_lib.h>


// C++ includes

namespace fest {

template <int dim>
class Quadrature
{

public:
    Quadrature(std::shared_ptr< dealii::Quadrature<dim> > quad_space,
               std::shared_ptr< dealii::Quadrature<1> > quad_time);
    
    std::shared_ptr< dealii::Quadrature<dim> > spatial();
    std::shared_ptr< dealii::Quadrature<1> > temporal();

private:
    std::shared_ptr< dealii::Quadrature<dim> > _quad_space;
    std::shared_ptr< dealii::Quadrature<1> > _quad_time;

};




}

#endif