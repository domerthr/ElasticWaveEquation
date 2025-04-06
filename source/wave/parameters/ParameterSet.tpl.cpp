// PROJECT includes
#include <wave/parameters/ParameterSet.hh>

// DEALI.II includes
#include <deal.II/base/exceptions.h>


// C++ includes


namespace wave {

ParameterSet::ParameterSet(std::shared_ptr< dealii::ParameterHandler > handler) {
    Assert(handler.use_count(), dealii::ExcNotInitialized());

    junker = handler->get_bool("Junker Boundary");
    dim = static_cast<unsigned int>(handler->get_integer("dim"));
    example = static_cast<unsigned int>(handler->get_integer("example"));
    s = static_cast<unsigned int>(handler->get_integer("s"));
    r = static_cast<unsigned int>(handler->get_integer("r"));


    handler->enter_subsection("Space mesh Specification");
    {
        TriaClass = handler->get("TriaClass");
        TriaOptions = handler->get("TriaOptions");

        global_refinement_space = static_cast <unsigned int>(
                                      handler->get_integer("global refinement space"));
        refine_space = handler->get_bool("refine space");
    }
    handler->leave_subsection();

    handler->enter_subsection("Time mesh Specification");
    {
        t0 = handler->get_double("initial time");
        T = handler->get_double("final time");
        M = handler->get_integer("time step");
        global_refinement_time = static_cast <unsigned int>(
                                     handler->get_integer("global refinement time"));
        refine_time = handler->get_bool("refine time");
        max_refinements_time = static_cast <unsigned int>(
                                   handler->get_integer("max refinements time"));

    }
    handler->leave_subsection();

    handler->enter_subsection("Parameter Specification");
    {
        lame_coefficient_mu = handler->get_double("mu");
        poisson_ratio_nu = handler->get_double("nu");
        rho = handler->get_double("rho");
        omega_dot = handler->get_double("omega_dot");

    }
    handler->leave_subsection();
}
}