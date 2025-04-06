// PROJECT includes
#include <wave/FinalValue/FinalValue_Example1.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

template <int dim>
double FinalValueExample1<dim>::value(const dealii::Point<dim> &p, 
                                      const unsigned int component) const {
    
    
    if (component == 0)
        return std::sin(_T) * p[0] * (1. - p[0]);
    else if (component == 1)
        return std::cos(_T) * p[0] * (1. - p[0]);
    else
        return -1.;    
    

}

template <int dim>
void FinalValueExample1<dim>::vector_value(const dealii::Point<dim> &p,
                                             dealii::Vector<double> &values) const {
  for (unsigned int c = 0; c < 2; ++c)
    values(c) = FinalValueExample1<dim>::value(p, c);
}
} // namespace wave 

// Implementierung für die Dimension 1
#include "FinalValue_Example1.inst.in"