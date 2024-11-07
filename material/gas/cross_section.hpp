#ifndef flick_material_cross_section
#define flick_material_cross_section

#include "../../numeric/function.hpp"

namespace flick {
  class basic_cross_section {
  protected:
    const std::string data_path = "material/gas/cross_section_input/";
  };
  
  class no2_cross_section : public basic_cross_section {
     pp_function c = read<pp_function>(data_path+"cross_section_no2.txt");
  public:
    double value(double wavelength) {
      return c.value(wavelength);
    }
    double longest() {
      return c.x().back();
    }
  };
  
  class o3_cross_section : public basic_cross_section {
    pp_function c0 = read<pp_function>(data_path+"ozone_cross_section_293K.txt");
    pl_function s = read<pl_function>(data_path+"ozone_temperature_dependence.txt");
    double T0 = 293;
  public:
    double value(double wavelength, double temperature) {
      return c0.value(wavelength)*exp(s.value(wavelength)*(temperature-T0));
    }
    double longest() {
      return c0.x().back();
    }
  };
}

#endif
