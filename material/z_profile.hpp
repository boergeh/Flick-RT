#ifndef flick_material_z_profile
#define flick_material_z_profile

#include "material.hpp"

namespace flick {
namespace material { 
  class z_profile : public base {
  protected:
    iop_z_profile a_profile_;
    iop_z_profile s_profile_;
  public:
    double absorption_coefficient() {
      return a_profile_.value(pose().position().z());
    }
    double scattering_coefficient() {
      return s_profile_.value(pose().position().z());
    }
    double absorption_optical_depth(double distance) {
      return a_profile_.optical_depth(pose(),distance);
    }
    double scattering_optical_depth(double distance) {
      return s_profile_.optical_depth(pose(),distance);
    }
    double absorption_distance(double absorption_optical_depth) {
      return s_profile_.distance(pose(),absorption_optical_depth);
    }
    double scattering_distance(double scattering_optical_depth) {
      return s_profile_.distance(pose(),scattering_optical_depth);
    }
  };
}
}

#endif
