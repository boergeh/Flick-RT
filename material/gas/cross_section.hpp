#ifndef flick_material_cross_section
#define flick_material_cross_section

#include "../../numeric/function.hpp"

namespace flick {
  class basic_cross_section {
  protected:
    std::vector<pp_function> cross_sections_;
    const std::string data_path = "material/gas/cross_section_input/";
  public:
    double longest() {
      return cross_sections_[0].x().back();
    }
  protected:
    void add_spectrum(const std::string& file_name) {
      cross_sections_.emplace_back(read<pp_function>(data_path+file_name));
    }
  };
  
  class no2_cross_section : public basic_cross_section {
  public:
    no2_cross_section() {
      add_spectrum("cross_section_no2.txt");
    }
    double value(double wavelength) {
      return cross_sections_[0].value(wavelength);
    }
  };
  
  class o3_cross_section : public basic_cross_section {
    std::vector<int> temperatures_{221,241,273};
  public:
    o3_cross_section() {
      for (size_t i=0; i<temperatures_.size(); ++i) {
	add_spectrum("cross_section_o3_"+ std::to_string(temperatures_[i])+"K.txt");
      }
    }
    double value(double wavelength, double temperature) {
      pp_function f;
      for (size_t i=0; i<cross_sections_.size(); ++i) {
	f.append({static_cast<double>(temperatures_[i]), cross_sections_[i].value(wavelength)});
      }
      return f.value(temperature); 
    }
  };
}

#endif
