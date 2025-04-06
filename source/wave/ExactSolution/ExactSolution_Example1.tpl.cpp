// PROJECT includes
#include <wave/ExactSolution/ExactSolution_Example1.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {

template <int dim>
double ExactSolutionExample1<dim>::value(const dealii::Point<dim> &p , 
                                         const unsigned int component) const {
    
    const double t = this->get_time();

    if (component == 0)
        return std::sin(t) * p[0] * (1 - p[0]); // u(t,x) = sin(t)x(1 - x)
    else
        return std::cos(t) * p[0] * (1 - p[0]); // v(t,x) = cos(t)x(1 - x)
        

}

template <int dim>
void ExactSolutionExample1<dim>::vector_value(const dealii::Point<dim> &p,
                                              dealii::Vector<double> &values) const
{
    for (unsigned int c = 0; c < 2; ++c)
        values(c) = ExactSolutionExample1<dim>::value(p,c);

}

}

#include "ExactSolution_Example1.inst.in"