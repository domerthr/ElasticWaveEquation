#ifndef __InitialValue_Example3_tpl_hh
#define __InitialValue_Example3_tpl_hh

// PROJECT includes

// DEAL.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>


// C++ includes

namespace wave {

template <int dim>
class InitialValueExample3 : public dealii::Function<dim>
{


public:
    InitialValueExample3() : dealii::Function<dim> (2 * dim) {};

    virtual ~InitialValueExample3() = default;
    
    virtual double value(const dealii::Point<dim> &p, 
                         const unsigned int component) const override;

    virtual void vector_value(const dealii::Point<dim> &p, 
                              dealii::Vector<double> &value) const override;
    

};
}

#endif