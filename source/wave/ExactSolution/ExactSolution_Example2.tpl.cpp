// PROJECT includes
#include <wave/ExactSolution/ExactSolution_Example2.tpl.hh>

// DEAL.II includes

// iDEAL.II includes


// C++ includes

namespace wave {

template <int dim>
double ExactSolutionExample2<dim>::value(const dealii::Point<dim> &p , 
                                         const unsigned int component) const {
    
    const double t = this->get_time();
    const double x = p[0];

    if (component == 0)
        return std::sin(M_PI * t) * std::sin(M_PI * x);  // u(t,x) = sin(πt) sin(πx)
    else
        return M_PI * std::cos(M_PI * t) * std::sin(M_PI * x);  // v(t,x) = π cos(πt) sin(πx)
        

}

template <int dim>
void ExactSolutionExample2<dim>::vector_value(const dealii::Point<dim> &p,
                                              dealii::Vector<double> &values) const
{
    for (unsigned int c = 0; c < 2; ++c)
        values(c) = ExactSolutionExample2<dim>::value(p,c);

}

}

#include "ExactSolution_Example2.inst.in"