// PROJECT includes
#include <wave/FinalValue/FinalValue_Example4.tpl.hh>

// DEAL.II includes

// iDEAL.II includes


// C++ includes

namespace wave {
    

template <int dim>
double FinalValueExample4<dim>::value(const dealii::Point<dim> &p, 
                                      const unsigned int component) const {
    
    const double x = p[0];
    if (component == 0)
        return std::sin(_T) * std::sin(M_PI * x);
    else if (component == 1)
        return std::cos(_T) * std::sin(M_PI * x);
    else
        return -1.;    
    

}

template <int dim>
void FinalValueExample4<dim>::vector_value(const dealii::Point<dim> &p,
                                             dealii::Vector<double> &values) const {
  for (unsigned int c = 0; c < 2; ++c)
    values(c) = FinalValueExample4<dim>::value(p, c);
}
} // namespace wave 

#include "FinalValue_Example4.inst.in"