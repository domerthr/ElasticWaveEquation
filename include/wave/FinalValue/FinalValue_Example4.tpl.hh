#ifndef __FinalValue_Example4_tpl_hh
#define __FinalValue_Example4_tpl_hh

// PROJECT includes

// DEAL.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>


// C++ includes


namespace wave {

template <int dim>
class FinalValueExample4 : public dealii::Function<dim>
{


public:
    FinalValueExample4(double T) : dealii::Function<dim> (2 * dim) {_T = T;};

    virtual ~FinalValueExample4() = default;
    
    virtual double value(const dealii::Point<dim> &p, 
                         const unsigned int component) const override;

    virtual void vector_value(const dealii::Point<dim> &p, 
                              dealii::Vector<double> &value) const override;

private:
    double _T;
    

};

}

#endif