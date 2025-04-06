// PROJECT includes
#include <wave/InitialValue/InitialValue_Example1.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

template <int dim>
double InitialValueExample1<dim>::value(const dealii::Point<dim> &p, 
                                        const unsigned int component) const {
    
    
    if (component == 0)
        return 0.;
    else if (component == 1)
        return p[0] * (1. - p[0]);
    else
        return -1.;    
    

}

template <int dim>
void InitialValueExample1<dim>::vector_value(const dealii::Point<dim> &p,
                                             dealii::Vector<double> &values) const {
    for (unsigned int c = 0; c < 2; ++c)
        values(c) = InitialValueExample1<dim>::value(p, c);
}

} // namespace wave 


#include "InitialValues_Example1.inst.in"