// PROJECT includes
#include <fest/base/refinement_table.hh>

// DEAL.II includes

// C++ includes


namespace fest {

RefinementTable::RefinementTable() {
    _max_refinement_cycles = 0;

}
void RefinementTable::print_row() {
    const int bar_width = 26;
    float progress = static_cast<float>(_current) / _total;
    int pos = static_cast<int>(bar_width * progress);
    std::ostringstream out;
    out << std::fixed << std::setprecision(4) << _time; // Genau 5 Nachkommastellen
    std::string result = out.str();
    
    std::cout << "\r" << std::string(7 - std::to_string(_refinement_cycle).length(), ' ')
              << _refinement_cycle 
              << std::string(5, ' ') << "|"
              << std::string(6 - std::to_string(_dofs_space).length(), ' ')
              << _dofs_space
              << " |" << std::string(5 - std::to_string(_dofs_time).length(), ' ')
              << _dofs_time
              << " |";
    std::cout << " [";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "]"<< std::string(4 - std::to_string(int(progress * 100.0)).length(), ' ') << int(progress * 100.0) << "% " << "|";
    if (_time == 1e12) std::cout << std::flush;
    else {
        
        std::cout << std::string(10-result.length(), ' ')
                  << result << "s" << "\n";
        _time = 1e12;
        _current = 0;
    }
}

void RefinementTable::print_header() {
    std::cout << std::string(76, '-') << "\n";
    std::cout << std::string(12, ' ') << "|"
              << std::string(5, ' ') << "dofs" << std::string(5, ' ') << "|"
              << std::string(34, ' ') << " |"
              << "\n";
    std::cout << " Refinement " << "|"
              << " space " << "|" << " time " <<  "|"
              << std::string(12, ' ') << " Progress " << std::string(12, ' ') << " |"
              << std::string(4, ' ') <<"Time" << "\n";
    std::cout << std::string(76, '-') << "\n";
}

template <typename T>
void RefinementTable::add_value(const std::string &column, const T value) {
    if (column.compare("dofs (space)") == 0) _dofs_space = value;
    else if (column.compare("total") == 0) _total = value;
    else if (column.compare("dofs (time)") == 0) _dofs_time = value;
    else if (column.compare("time") == 0) _time = value;
    else if (column.compare("current slab") == 0) _current = value;
    else if (column.compare("refinement cycle") == 0) _refinement_cycle = value;

    
}

}

#include "refinement_table.inst.in"
