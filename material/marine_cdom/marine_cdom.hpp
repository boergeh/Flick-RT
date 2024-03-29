#ifndef flick_material_marine_cdom
#define flick_material_marine_cdom

#include "../material.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class marine_cdom : public base {
    pl_function a_;
    double a440_;
  public:
    marine_cdom(const std::string& name, double scaling_factor) {
      const std::string path = "/material/marine_cdom/iop_tables";
      a_ = read<pl_function>(add_path_if_exists(name+"_a.txt",{"./",path}));
      a_.scale_y(scaling_factor);
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
