#ifndef __TriaGenerator_tpl_hh
#define __TriaGenerator_tpl_hh

// PROJECT includes

// DEALI.II includes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/point.h>

// C++ includes
#include <fstream>
#include <iostream>
#include <string>

namespace wave {
namespace grid {


template <int dim>
class TriaGenerator {

public:
    TriaGenerator() = default;
    void generate(std::shared_ptr< dealii::Triangulation<dim> > tria,
                  const unsigned int example);

};

   
}}

#endif