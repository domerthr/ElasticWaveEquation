#ifndef __fe_tpl_hh
#define __fe_tpl_hh

// PROJECT includes


// DEAL.II includes
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>


// C++ includes

namespace fest {
namespace cgcg {

template <int dim>
class FiniteElement {

public:

    FiniteElement(std::shared_ptr< dealii::FiniteElement<dim>> fe_space, 
                  const unsigned int r);

    const unsigned int dofs_per_cell;

    std::shared_ptr< dealii::FiniteElement<dim> >   spatial();
    std::shared_ptr< dealii::FiniteElement<1> >     temporal();

private:
    std::shared_ptr< dealii::FiniteElement<dim> >   _fe_space;
    std::shared_ptr< dealii::FiniteElement<1> >     _fe_time;


};
} // namespace cgcg


namespace cg {
    
    
template <int dim>
class FiniteElement
{
public:

    enum support_type
    {
        /**Support points based on QGauss<1>*/
        Legendre,
        /**for dG(r), r>0: Support points based on QGaussLobatto<1>*/
        Lobatto,
        /**Support points based on left QGaussRadau<1>*/
        RadauLeft,
        /**Support points based on right QGaussRadau<1>*/
        RadauRight
    };


    FiniteElement(std::shared_ptr<dealii::FiniteElement<dim>> fe_space,
                  const unsigned int                          r,
                  support_type type = support_type::Lobatto);

    std::shared_ptr<dealii::FiniteElement<dim>>
    spatial();

    std::shared_ptr<dealii::FiniteElement<1>>
    temporal_trial();

    std::shared_ptr<dealii::FiniteElement<1>>
    temporal_test();

    const unsigned int dofs_per_cell_trial;
    const unsigned int dofs_per_cell_test;

    support_type
    type();


private:
    std::shared_ptr<dealii::FiniteElement<dim>> _fe_space;
    std::shared_ptr<dealii::FiniteElement<1>>   _fe_time_trial;
    std::shared_ptr<dealii::FiniteElement<1>>   _fe_time_test;
    support_type                                _type;

};



} // namespace cg
} // namespace fest

#endif