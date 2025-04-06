// PROJECT includes
#include <wave/FinalValue/FinalValue_Selector.tpl.hh>
#include <wave/FinalValue/FinalValues.hh>

namespace wave {
namespace final_value {
    
template <int dim>
void Selector<dim>::create_function(const unsigned int example,
                                    std::shared_ptr< dealii::Function<dim>> &function,
                                    double T) const {

    if (example == 1) {

        function = std::make_shared< FinalValueExample1<dim> >(T);
        return;
    }
    else if (example == 2){

        function = std::make_shared< FinalValueExample2<dim> >(T);
        return;
    }
    else if (example == 4){

        function = std::make_shared< FinalValueExample4<dim> >(T);
        return;
    }
    else {

        function = nullptr;
        return;
    }
    

    AssertThrow(false,
		        dealii::ExcMessage("final_value function unknown, please check your input file data.")
	);
    
}

} // namespace initial_value
}

#include "FinalValue_Selector.inst.in"