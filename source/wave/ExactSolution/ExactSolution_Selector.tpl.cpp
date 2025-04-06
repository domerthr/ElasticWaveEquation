// PROJECT includes
#include <wave/ExactSolution/ExactSolution_Selector.tpl.hh>
#include <wave/ExactSolution/ExactSolutions.hh>

// DEAL.II includes



// C++ includes

namespace wave {
namespace exact_solution {
    
template <int dim>
void Selector<dim>::create_function(const unsigned int example,
                                    std::shared_ptr< dealii::Function<dim>> &function) const {

    if (example == 1) {

        function = std::make_shared< ExactSolutionExample1<dim> >();
        return;
    }
    else if (example == 2) {

        function = std::make_shared< ExactSolutionExample2<dim> >();
        return;
    }
    else if (example == 4) {

        function = std::make_shared< ExactSolutionExample4<dim> >();
        return;
    }
    else {

        function = nullptr;
        return;
    }

    AssertThrow(false,
		        dealii::ExcMessage("exact solution function unknown, please check your input file data.")
	);



    
}

} // namespace exact_solution
}

#include "ExactSolution_Selector.inst.in"