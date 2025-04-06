// PROJECT includes
#include <wave/ExactSolution/ExactSolution_Example4.tpl.hh>

// DEAL.II includes

// iDEAL.II includes


// C++ includes

namespace wave {

template <int dim>
double ExactSolutionExample4<dim>::value(const dealii::Point<dim> &p , 
                                         const unsigned int component) const {
    
    const double t = this->get_time();
    const double x = p[0];

    if (component == 0)
        return std::sin(t) * std::sin(M_PI * x);  // u(t,x) = sin(t) sin(πx)
    else
        return std::cos(t) * std::sin(M_PI * x);  // v(t,x) = cos(t) sin(πx)
        

}

template <int dim>
void ExactSolutionExample4<dim>::vector_value(const dealii::Point<dim> &p,
                                              dealii::Vector<double> &values) const
{
    for (unsigned int c = 0; c < 2; ++c)
        values(c) = ExactSolutionExample4<dim>::value(p,c);

}

}

#include "ExactSolution_Example4.inst.in"