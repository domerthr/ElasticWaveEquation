// PROJECT includes
#include <wave/FinalValue/FinalValue_Example2.tpl.hh>

// DEAL.II includes

// iDEAL.II includes


// C++ includes

namespace wave {
    

template <int dim>
double FinalValueExample2<dim>::value(const dealii::Point<dim> &p, 
                                      const unsigned int component) const {
    
    const double x = p[0];
    if (component == 0)
        return std::sin(M_PI * _T) * std::sin(M_PI * x);
    else if (component == 1)
        return M_PI * std::cos(M_PI * _T) * std::sin(M_PI * x);
    else
        return -1.;    
    

}

template <int dim>
void FinalValueExample2<dim>::vector_value(const dealii::Point<dim> &p,
                                             dealii::Vector<double> &values) const {
  for (unsigned int c = 0; c < 2; ++c)
    values(c) = FinalValueExample2<dim>::value(p, c);
}
} // namespace wave 

#include "FinalValue_Example2.inst.in"