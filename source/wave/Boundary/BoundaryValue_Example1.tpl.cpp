// PROJECT includes
#include <wave/Boundary/BoundaryValue_Example1.tpl.hh>


// DEAL.II includes



// C++ includes

namespace wave {

template <int dim>
double BoundaryValueExample1<dim>::value([[maybe_unused]] const dealii::Point<dim> &p, 
                                         [[maybe_unused]] const unsigned int component) const {
    
    Assert(component < this->n_components, dealii::ExcIndexRange(component, 0, this->n_components));
    
    return 0.;

}

template <int dim>
void BoundaryValueExample1<dim>::vector_value(const dealii::Point<dim> &p,
                                              dealii::Vector<double> &values) const {
    for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = BoundaryValueExample1<dim>::value(p, c);
}

} // namespace


#include "BoundaryValue_Example1.inst.in"