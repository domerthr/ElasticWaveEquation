// PROJECT includes
#include <wave/Boundary/BoundaryValue_Selector.tpl.hh>
#include <wave/Boundary/BoundaryValues.hh>


// DEAL.II includes


// C++ includes
#include <iostream>

namespace wave {
namespace boundary_value {
    
template <int dim>
void Selector<dim>::create_function(const unsigned int example,
                                    std::shared_ptr< dealii::Function<dim>> &function) const {

    if (example == 1) {


        function = std::make_shared< BoundaryValueExample1<dim> >();
        return;
    }
    else if (example == 2) {

        function = std::make_shared< BoundaryValueExample2<dim> >();
        return;
    }
    else if (example == 3) {

        function = std::make_shared< BoundaryValueExample3<dim> >();
        return;
    }
    else if (example == 4) {

        function = std::make_shared< BoundaryValueExample4<dim> >();
        return;
    }
    

    AssertThrow(false,
		        dealii::ExcMessage("boundary_value function unknown, please check your input file data.")
	);
    
}

} // namespace boundary_value
} // namespace wave

#include "BoundaryValue_Selector.inst.in"