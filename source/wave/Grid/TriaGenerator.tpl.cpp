// PROJECT includes
#include <wave/Grid/TriaGenerator.tpl.hh>

// DEALI.II includes

// C++ includes

namespace wave {
namespace grid {


template <int dim>
void
TriaGenerator<dim>::generate(std::shared_ptr< dealii::Triangulation<dim> > tria,
                             const unsigned int example) {
    

    if (example == 1) {
        dealii::GridGenerator::hyper_rectangle(*tria, dealii::Point<dim>(0.), dealii::Point<dim>(1.));


        auto cell(tria->begin_active());
        auto endc(tria->end());

        // setting boundary id=0
        for (; cell != endc; ++cell) {
            for (unsigned int face(0); face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
                if (cell->face(face)->at_boundary()) {
                    cell->face(face)->set_boundary_id(0);    
                }
            }
        }
    }
    else if (example == 2) {
        dealii::GridGenerator::hyper_rectangle(*tria, dealii::Point<dim>(0.), dealii::Point<dim>(2.));


        auto cell(tria->begin_active());
        auto endc(tria->end());

        // setting boundary id=0
        for (; cell != endc; ++cell) {
            for (unsigned int face(0); face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
                if (cell->face(face)->at_boundary()) {
                    cell->face(face)->set_boundary_id(0);    
                }
            }
        }
    }

    else if (example == 3) {
        dealii::Point<dim> bottom_left, upper_right;
        upper_right[0] = 5.;
        upper_right[1] = .5;
        dealii::GridGenerator::hyper_rectangle(*tria, bottom_left, upper_right);
        
        for (auto &face : tria->active_face_iterators())
            if (face->at_boundary())
                if (std::fabs(face->center()(1) - (.5)) < 1e-12 ||
                    std::fabs(face->center()(1) - (0.)) < 1e-12 ||
                    std::fabs(face->center()(0) - (5.)) < 1e-12)
                face->set_boundary_id(1);
    

    }
    else if (example == 4) {
        dealii::GridGenerator::hyper_rectangle(*tria, dealii::Point<dim>(0.), dealii::Point<dim>(1.));


        auto cell(tria->begin_active());
        auto endc(tria->end());

        // setting boundary id=0
        for (; cell != endc; ++cell) {
            for (unsigned int face(0); face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
                if (cell->face(face)->at_boundary()) {
                    cell->face(face)->set_boundary_id(0);    
                }
            }
        }
    }
        

    

}

   
}}

#include "TriaGenerator.inst.in"
