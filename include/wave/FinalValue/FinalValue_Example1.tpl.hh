#ifndef __FinalValue_Example1_tpl_hh
#define __FinalValue_Example1_tpl_hh

// PROJECT includes

// DEAL.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>


// C++ includes


namespace wave {

template <int dim>
class FinalValueExample1 : public dealii::Function<dim>
{


public:
    FinalValueExample1(double T) : dealii::Function<dim> (2 * dim) {_T = T;};

    virtual ~FinalValueExample1() = default;
    
    virtual double value(const dealii::Point<dim> &p, 
                         const unsigned int component) const override;

    virtual void vector_value(const dealii::Point<dim> &p, 
                              dealii::Vector<double> &value) const override;
private:
    double _T;

};

}

#endif