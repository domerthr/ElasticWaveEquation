#ifndef __SpaceTime_dof_tools_hh
#define __SpaceTime_dof_tools_hh

// PROJECT includes
#include <fest/dofs/SpaceTime_dof_handler.tpl.hh>
#include <fest/fe/fe.tpl.hh>


// DEAL.II includes
#include <deal.II/base/types.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>


// C++ includes



namespace fest {
namespace cgcg {
namespace DoFTools {

template <int dim>
void make_sparsity_pattern(fest::cgcg::DoFHandler<dim> &dof, dealii::DynamicSparsityPattern &dsp,
                           std::shared_ptr< dealii::AffineConstraints<double> > space_constraints = std::make_shared<dealii::AffineConstraints<double>>(),
                           const bool keep_constrained_dofs = true) {

    
    // spartial sparsity pattern
    dealii::DynamicSparsityPattern space_dsp(dof.n_dofs_space());
    dealii::DoFTools::make_sparsity_pattern(*dof.spatial(), space_dsp, *space_constraints, keep_constrained_dofs);

    ////dealii::SparsityPattern space_sparsity_pattern;
    ////space_sparsity_pattern.copy_from(space_dsp);
    ////std::ofstream out_space_sparsity("../space_sparsity_pattern.svg");
    ////space_sparsity_pattern.print_svg(out_space_sparsity);


    // temporal sparsity pattern
    dealii::DynamicSparsityPattern time_dsp(dof.n_dofs_time());

    dealii::DoFTools::make_sparsity_pattern(*dof.temporal(), time_dsp);
            
    for (auto &space_entry : space_dsp) {
        for (auto &time_entry : time_dsp) {
            dsp.add(time_entry.row() * dof.n_dofs_space() + space_entry.row(),         // test function
                    time_entry.column() * dof.n_dofs_space() + space_entry.column()    // trial function
                    );
        }
    }
}


} // namespace DoFTools
} // namespace cgcg

namespace cg {
namespace DoFTools {
// namespace internal
// {

// template <int dim>
// void
// upwind_temporal_pattern(fest::cg::slab::DoFHandler<dim> &dof,
//                         dealii::DynamicSparsityPattern  &time_dsp)
// {
//     dealii::DoFTools::make_sparsity_pattern(*dof.temporal(), time_dsp);
//     if (dof.fe_support_type() == fest::cg::FiniteElement<dim>::support_type::Lobatto)
//     {
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time(); ii < dof.n_dofs_time(); ii += dof.dofs_per_cell_time())
//         {
//             time_dsp.add(ii, ii - 1);
//         }
//     }
//     else if (dof.fe_support_type() == fest::cg::FiniteElement<dim>::support_type::RadauLeft)
//     {
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time();
//             ii < dof.n_dofs_time();
//             ii += dof.dofs_per_cell_time())
//         {
//             for (dealii::types::global_dof_index l = 0;
//                  l < dof.dofs_per_cell_time();
//                  l++)
//             {
//                 time_dsp.add(ii, ii - l - 1);
//             }
//         }
//     }
//     else if (dof.fe_support_type() == fest::cg::FiniteElement<dim>::support_type::RadauRight)
//     {
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time();
//              ii < dof.n_dofs_time();
//              ii += dof.dofs_per_cell_time())
//         {
//             for (dealii::types::global_dof_index k = 0;
//                  k < dof.dofs_per_cell_time();
//                  k++)
//             {
//                 time_dsp.add(ii + k, ii - 1);
//             }
//         }
//     }
//     else
//     {
//         // go over first DoF of each cell
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time();
//              ii < dof.n_dofs_time();
//              ii += dof.dofs_per_cell_time())
//         {
//             // row offset
//             for (dealii::types::global_dof_index k = 0;
//                  k < dof.dofs_per_cell_time();
//                  k++)
//             {
//                 for (dealii::types::global_dof_index l = 0;
//                      l < dof.dofs_per_cell_time();
//                      l++)
//                 {
//                     time_dsp.add(ii + k, ii - l - 1);
//                 }
//             }
//         }
//     }
// }

// template <int dim>
// void
// downwind_temporal_pattern(fest::cg::slab::DoFHandler<dim> &dof,
//                           dealii::DynamicSparsityPattern  &time_dsp)
// {
//     dealii::DoFTools::make_sparsity_pattern(*dof.temporal(), time_dsp);
//     if (dof.fe_support_type() == fest::cg::FiniteElement<dim>::support_type::Lobatto)
//     {
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time();
//              ii < dof.n_dofs_time();
//              ii += dof.dofs_per_cell_time())
//         {
//             time_dsp.add(ii - 1, ii);
//         }
//     }
//     else if (dof.fe_support_type() == fest::cg::FiniteElement<dim>::support_type::RadauRight)
//     {
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time();
//              ii < dof.n_dofs_time();
//              ii += dof.dofs_per_cell_time())
//         {
//             for (dealii::types::global_dof_index k = 0;
//                     k < dof.dofs_per_cell_time();
//                     k++)
//             {
//                 time_dsp.add(ii - 1, ii + k);
//             }
//         }
//     }
//     else if (dof.fe_support_type() == fest::cg::FiniteElement<dim>::support_type::RadauLeft)
//     {
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time();
//              ii < dof.n_dofs_time();
//              ii += dof.dofs_per_cell_time())
//         {
//             for (dealii::types::global_dof_index l = 0;
//                  l < dof.dofs_per_cell_time();
//                  l++)
//             {
//                 time_dsp.add(ii - l - 1, ii);
//             }
//         }
//     }
//     else
//     {
//         // go over first DoF of each cell
//         for (dealii::types::global_dof_index ii = dof.dofs_per_cell_time();
//              ii < dof.n_dofs_time();
//              ii += dof.dofs_per_cell_time())
//         {
//             // row offset
//             for (dealii::types::global_dof_index k = 0;
//                  k < dof.dofs_per_cell_time();
//                  k++)
//             {
//                 for (dealii::types::global_dof_index l = 0;
//                      l < dof.dofs_per_cell_time();
//                      l++)
//                 {
//                     time_dsp.add(ii - l - 1, ii + k);
//                 }
//             }
//         }
//     }
// }
// } // namespace internal
    
// template <int dim>
// void
// make_upwind_sparsity_pattern(fest::cg::slab::DoFHandler<dim> &dof,
//                              dealii::DynamicSparsityPattern  &st_dsp,
//                              std::shared_ptr<dealii::AffineConstraints<double>> space_constraints =
//                                     std::make_shared<dealii::AffineConstraints<double>>(),
//                              const bool keep_constrained_dofs = true)
// {
//     dealii::DynamicSparsityPattern space_dsp(dof.n_dofs_space());
//     dealii::DoFTools::make_sparsity_pattern(*dof.spatial(),
//                                             space_dsp,
//                                             *space_constraints,
//                                             keep_constrained_dofs);

//     dealii::DynamicSparsityPattern time_dsp(dof.n_dofs_time());
//     internal::upwind_temporal_pattern(dof, time_dsp);

//     for (auto &space_entry : space_dsp) {
//         for (auto &time_entry : time_dsp) {
//             st_dsp.add(time_entry.row() * dof.n_dofs_space() + space_entry.row(), // test function
//                        time_entry.column() * dof.n_dofs_space() + space_entry.column() // trial function
//             );
//         }
//     }
// }

// template <int dim>
// void
// make_downwind_sparsity_pattern(fest::cg::slab::DoFHandler<dim> &dof,
//                                dealii::DynamicSparsityPattern  &st_dsp,
//                                std::shared_ptr<dealii::AffineConstraints<double>> space_constraints =
//                                     std::make_shared<dealii::AffineConstraints<double>>(),
//                                const bool keep_constrained_dofs = true) {


//     dealii::DynamicSparsityPattern space_dsp(dof.n_dofs_space());
//     dealii::DoFTools::make_sparsity_pattern(*dof.spatial(),
//                                             space_dsp,
//                                             *space_constraints,
//                                             keep_constrained_dofs);

//     dealii::DynamicSparsityPattern time_dsp(dof.n_dofs_time());
//     internal::downwind_temporal_pattern(dof, time_dsp);

//     for (auto &space_entry : space_dsp) {
//         for (auto &time_entry : time_dsp) {
//             st_dsp.add(time_entry.row() * dof.n_dofs_space() + space_entry.row(), // test function
//                        time_entry.column() * dof.n_dofs_space() + space_entry.column() // trial function
//             );
//         }
//     }
// }

template <int dim>
void make_sparsity_pattern(fest::cg::slab::DoFHandler<dim> &dof, dealii::DynamicSparsityPattern &dsp,
                           std::shared_ptr< dealii::AffineConstraints<double> > space_constraints = std::make_shared<dealii::AffineConstraints<double>>(),
                           const bool keep_constrained_dofs = true) {

    
    // spartial sparsity pattern
    dealii::DynamicSparsityPattern space_dsp(dof.n_dofs_space());
    dealii::DoFTools::make_sparsity_pattern(*dof.spatial(), space_dsp, *space_constraints, keep_constrained_dofs);

    ////dealii::SparsityPattern space_sparsity_pattern;
    ////space_sparsity_pattern.copy_from(space_dsp);
    ////std::ofstream out_space_sparsity("../space_sparsity_pattern.svg");
    ////space_sparsity_pattern.print_svg(out_space_sparsity);


    // temporal sparsity pattern
    dealii::DynamicSparsityPattern time_dsp(dof.n_dofs_time_trial());

    dealii::DoFTools::make_sparsity_pattern(*dof.temporal_trial(), time_dsp);
            
    for (auto &space_entry : space_dsp) {
        for (auto &time_entry : time_dsp) {
            dsp.add(time_entry.row() * dof.n_dofs_space() + space_entry.row(),         // test function
                    time_entry.column() * dof.n_dofs_space() + space_entry.column()    // trial function
                    );
        }
    }
}

} // namespace DoFTools
} // namespace cg


}

#endif