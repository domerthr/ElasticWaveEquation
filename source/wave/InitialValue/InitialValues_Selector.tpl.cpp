// PROJECT includes
#include <wave/InitialValue/InitialValue_Selector.tpl.hh>
#include <wave/InitialValue/InitialValues.hh>

// DEAL.II includes



// C++ includes

namespace wave {
    

namespace initial_value {
    
template <int dim>
void Selector<dim>::create_function(const unsigned int example,
                                    std::shared_ptr< dealii::Function<dim>> &function) const {

    if (example == 1) {

        function = std::make_shared< InitialValueExample1<dim> >();
        return;
    }
    else if (example == 2) {

        function = std::make_shared< InitialValueExample2<dim> >();
        return;
    }
    else if (example == 3) {

        function = std::make_shared< InitialValueExample3<dim> >();
        return;
    }
    else if (example == 4) {

        function = std::make_shared< InitialValueExample4<dim> >();
        return;
    }


    AssertThrow(false,
		        dealii::ExcMessage("initial_value function unknown, please check your input file data.")
	);
    
}

} // namespace initial_value
} // namespace wave


#include "InitialValues_Selector.inst.in"