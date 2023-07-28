#ifndef flick_material_cross_section
#define flick_material_cross_section

#include "../../numeric/function.hpp"

namespace flick {
  class o3_cross_section {
    std::vector<pp_function> cross_sections_;
    std::vector<int> temperatures_{221,241,273};
  public:
    o3_cross_section() {
      for (size_t i=0; i<temperatures_.size(); ++i) {
	std::string fname = "cross_section_o3_"+ std::to_string(temperatures_[i])+"K.txt";
	cross_sections_.emplace_back(read<pp_function>("material/gas/cross_section_input/"+fname));
      }
    }
    double value(double wavelength, double temperature) {
      pp_function f;
      for (size_t i=0; i<cross_sections_.size(); ++i) {
	f.append({static_cast<double>(temperatures_[i]), cross_sections_[i].value(wavelength)});
      }
      return f.value(temperature); 
    }
    double longest() {
      return cross_sections_[0].x().back();
    }
  };
}

#endif
