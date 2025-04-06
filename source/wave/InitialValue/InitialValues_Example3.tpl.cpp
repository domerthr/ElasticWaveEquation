// PROJECT includes
#include <wave/InitialValue/InitialValue_Example3.tpl.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

template <int dim>
double InitialValueExample3<dim>::value([[maybe_unused]] const dealii::Point<dim> &p, 
                                        [[maybe_unused]] const unsigned int component) const {
    
    return 0.;    
    

}

template <int dim>
void InitialValueExample3<dim>::vector_value(const dealii::Point<dim> &p,
                                             dealii::Vector<double> &values) const {
    for (unsigned int c = 0; c < this->n_components; ++c)
      values(c) = InitialValueExample3<dim>::value(p, c);
}

} // namespace wave 


#include "InitialValues_Example3.inst.in"