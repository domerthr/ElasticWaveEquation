#ifndef __SpaceTimeTria_tpl_hh
#define __SpaceTimeTria_tpl_hh

// PROJECT includes


// DEAL.II includes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>


// C++ includes
#include <list>
#include <memory>


namespace fest {
namespace cgcg {

template <int dim>
class Triangulation {

public:

    Triangulation();
    
    Triangulation(const Triangulation &other);
    
    std::shared_ptr< dealii::Triangulation<dim> > spatial();
    std::shared_ptr< dealii::Triangulation<1> > temporal();

    void generate(double startpoint, 
                  double endpoint);

    double startpoint();
    double endpoint();

    virtual void refine_global(const unsigned int n = 1, const unsigned int m = 1);

    unsigned int n_active_time_cells();
    unsigned int n_active_space_cells();
    unsigned int n_active_cells();


private:
    std::shared_ptr< dealii::Triangulation<1> >     _temporal_tria;
    std::shared_ptr< dealii::Triangulation<dim> >   _spatial_tria;
    double _startpoint;
    double _endpoint;

};

} // namespace cgcg
} // namespace fest

#endif