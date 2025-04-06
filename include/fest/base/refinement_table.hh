#ifndef __refinement_table_hh
#define __refinement_table_hh

// PROJECT includes

// DEAL.II includes

// C++ includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

namespace fest {

class RefinementTable {

public:
    RefinementTable();
    void print_header();
    void print_row();

    template <typename T>
    void add_value(const std::string &column, const T value);
    

private:
    unsigned int _max_refinement_cycles;
    unsigned int _refinement_cycle;
    unsigned int _current;
    unsigned int _total = 4;
    double _time = 1e12;
    unsigned int _dofs_time;
    unsigned int _dofs_space;


};
}

#endif