#ifndef __ExactSolution_Example4_tpl_hh
#define __ExactSolution_Example4_tpl_hh


// PROJECT includes

// DEAL.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>


// C++ includes

namespace wave {
template <int dim>
class ExactSolutionExample4 : public dealii::Function<dim>
{


public:
    ExactSolutionExample4() : dealii::Function<dim>(2 * dim) {};

    virtual ~ExactSolutionExample4() = default;

    virtual double value(const dealii::Point<dim> &p,
                         const unsigned int component = 0) const override;

    virtual void vector_value(const dealii::Point<dim> &p,
                              dealii::Vector<double> & value) const override;
    

};

}

#endif