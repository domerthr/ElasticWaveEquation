// PROJECT includes
#include <fest/fe/fe.tpl.hh>
#include "fest/base/quadrature_lib.tpl.hh"


// iDEAL.II includes
#include <ideal.II/base/quadrature_lib.hh>


// DEAL.II includes
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>


// C++ includes
#include <memory>

namespace fest {
namespace cgcg {

template <int dim>
FiniteElement<dim>::FiniteElement(std::shared_ptr< dealii::FiniteElement<dim> > fe_space,
                                  const unsigned int r)
    :   dofs_per_cell(fe_space->dofs_per_cell * (r + 1)) 
    ,   _fe_space(fe_space) 
{
    
    _fe_time = std::make_shared< dealii::FE_Q<1> > (r);
}

template <int dim>
std::shared_ptr< dealii::FiniteElement<dim> > FiniteElement<dim>::spatial() {

    return _fe_space;
}

template <int dim>
std::shared_ptr< dealii::FiniteElement<1> > FiniteElement<dim>::temporal() {
    return _fe_time;
}

} // namespace cgcg

namespace cg {

template <int dim>
FiniteElement<dim>::FiniteElement(std::shared_ptr<dealii::FiniteElement<dim>> fe_space,
                                  const unsigned int                          r,
                                  support_type                                type)
    :   dofs_per_cell_trial(fe_space->dofs_per_cell * (r + 1))
    ,   dofs_per_cell_test(fe_space->dofs_per_cell * r)
    ,   _fe_space(fe_space)
{
    _fe_time_trial = std::make_shared< dealii::FE_Q<1> > (r);

    if (type == support_type::Legendre ||
        (type == support_type::Lobatto && (r - 1) == 0)) {
        _fe_time_test = std::make_shared<dealii::FE_DGQArbitraryNodes<1>>(
          dealii::QGauss<1>(r));
    }
    else if (type == support_type::Lobatto) {
        _fe_time_test = std::make_shared<dealii::FE_DGQArbitraryNodes<1>>(
          dealii::QGaussLobatto<1>(r));
    }
    else if (type == support_type::RadauLeft) {
        _fe_time_test = std::make_shared<dealii::FE_DGQArbitraryNodes<1>>(
          idealii::QGaussRadau<1>(r, idealii::QGaussRadau<1>::left));
    }
    else {
        _fe_time_test = std::make_shared<dealii::FE_DGQArbitraryNodes<1>>(
          idealii::QGaussRadau<1>(r, idealii::QGaussRadau<1>::right));
    }
}

template <int dim>
std::shared_ptr<dealii::FiniteElement<dim>> FiniteElement<dim>::spatial() {
    return _fe_space;
}

template <int dim>
std::shared_ptr<dealii::FiniteElement<1>> FiniteElement<dim>::temporal_trial() {
    return _fe_time_trial;
}

template <int dim>
std::shared_ptr<dealii::FiniteElement<1>> FiniteElement<dim>::temporal_test() {
    return _fe_time_test;
}

template <int dim>
typename FiniteElement<dim>::support_type FiniteElement<dim>::type() {
    return _type;
}


} // namespace cg

} // namespace fest

#include "fe.inst.in"