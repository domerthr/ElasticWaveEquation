// PROJECT includes
#include <wave/parameters/ParameterHandler.hh>

// DEALI.II includes
#include <deal.II/base/parameter_handler.h>

// C++ includes

namespace wave {

ParameterHandler::ParameterHandler() {

    declare_entry(
        "Junker Boundary",
        "false",
        dealii::Patterns::Bool(),
        "Specification if Final Conditions should be used"
    );

    declare_entry(
        "dim",
        "1",
        dealii::Patterns::Integer(),
        "dim"
    );

    declare_entry(
        "example",
        "1",
        dealii::Patterns::Integer(),
        "example"
    );


    declare_entry(
        "s",
        "1",
        dealii::Patterns::Integer(),
        "s"
    );

    declare_entry(
        "r",
        "1",
        dealii::Patterns::Integer(),
        "r"
    );

    // The space mesh
    enter_subsection("Space mesh Specification");
    {
        declare_entry(
			"TriaClass",
			"invalid",
			dealii::Patterns::Anything()
		);
		
		declare_entry(
			"TriaOptions",
			"invalid",
			dealii::Patterns::Anything()
		);

        declare_entry(
			"global refinement space",
			"0",
			dealii::Patterns::Integer(),
			"Global refinements of the intial space mesh"
		);
        declare_entry(
            "refine space",
            "true",
            dealii::Patterns::Bool(),
            "Refine space mesh"
        );

    }
    leave_subsection();

    // The time mesh
    enter_subsection("Time mesh Specification");
    {
        declare_entry(
			"initial time",
			"0.",
			dealii::Patterns::Double(),
			"initial time t0"
		);
		
		declare_entry(
			"final time",
			"0.",
			dealii::Patterns::Double(),
			"final time T"
		);

        declare_entry(
			"time step",
			"1",
			dealii::Patterns::Integer(),
			"number of steps"
		);

        declare_entry(
			"global refinement time",
			"0",
			dealii::Patterns::Integer(),
			"Global refinements of the intial time mesh"
		);
        declare_entry(
            "refine time",
            "true",
            dealii::Patterns::Bool(),
            "Refine time mesh"
        );

        declare_entry(
			"max refinements time",
			"0",
			dealii::Patterns::Integer(),
			"Max refinements of the time mesh"
		);


    }
    leave_subsection();
    
    // Physical Parameters
    enter_subsection("Parameter Specification");
    {
        declare_entry(
            "mu",
            "1.0", 
            dealii::Patterns::Double(), 
            "Lamé Coefficient μ"
        );

        declare_entry(
            "nu",
            "1.0", 
            dealii::Patterns::Double(), 
            "Poisson ratio ν"
        );
  
        declare_entry(
            "rho", 
            "1.0", 
            dealii::Patterns::Double(), 
            "Density"
        );

        declare_entry(
            "omega_dot", 
            "1.0", 
            dealii::Patterns::Double(), 
            "Angular Acceleration"
        );

    }
    leave_subsection();

}}