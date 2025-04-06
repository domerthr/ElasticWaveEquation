// PROJECT includes
#include <wave/InitialValue/InitialValue_Example4.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

template <int dim>
double InitialValueExample4<dim>::value([[maybe_unused]] const dealii::Point<dim> &p, 
                                        const unsigned int component) const {
    
    
    if (component == 0)
        return 0.0;
    else if (component == 1)
        return std::sin(M_PI * p[0]);
    else
        return -1.;    
    

}

template <int dim>
void InitialValueExample4<dim>::vector_value(const dealii::Point<dim> &p,
                                             dealii::Vector<double> &values) const {
    for (unsigned int c = 0; c < 2; ++c)
        values(c) = InitialValueExample4<dim>::value(p, c);
}

} // namespace wave 


#include "InitialValues_Example4.inst.in"