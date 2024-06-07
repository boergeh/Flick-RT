#ifndef flick_material_marine_cdom
#define flick_material_marine_cdom

#include "../material.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class marine_cdom : public base {
    pl_function a_;
    double a440_;
    const double constant = 1;
  public:
    marine_cdom(const std::string& name, double scaling_factor) {
      std::string p = path()+"/material/marine_cdom/iop_tables";
      a_ = read<pl_function>(name+"_a.txt", p);
      a_.scale_y(scaling_factor);
      a_.add_extrapolation_points(constant);
    }
    double absorption_coefficient() const {
      const double to_nm = 1e9;
      return a_.value(wavelength()*to_nm);
    }
    double scattering_coefficient() const {
      return 0;
    }
  };
  
  template<int n>
  struct listable_marine_cdom : public marine_cdom {
    using marine_cdom::marine_cdom;
  };
}
}

#endif
