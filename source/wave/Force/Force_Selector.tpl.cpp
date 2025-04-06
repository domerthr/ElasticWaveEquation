// PROJECT includes
#include <wave/Force/Force_Selector.tpl.hh>
#include <wave/Force/Forces.hh>

// DEAL.II includes

// iDEAL.II includes

// C++ includes
#include <iostream>

namespace wave {
    

namespace force {
    
template <int dim>
void Selector<dim>::create_function(const unsigned int example,
                                    std::shared_ptr< dealii::TensorFunction<1, dim, double>> &function) const {

    if (example == 1) {

        
        function = std::make_shared< ForceExample1<dim> >();
        return;
    }
    else if (example == 2) {

        
        function = std::make_shared< ForceExample2<dim> >();
        return;
    }
    else if (example == 3) {

        
        function = std::make_shared< ForceExample3<dim> >();
        return;
    }
    else if (example == 4) {

        
        function = std::make_shared< ForceExample4<dim> >();
        return;
    }


    AssertThrow(false,
		        dealii::ExcMessage("force function unknown, please check your input file data.")
	);
    
}

} // namespace initial_value
} // namespace wave 


#include "Force_Selector.inst.in"