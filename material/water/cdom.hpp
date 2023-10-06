#ifndef flick_material_cdom
#define flick_material_cdom

#include "../material.hpp"

namespace flick {
namespace material {
  class cdom : public base {
    double abs_coef_440_;
    double slope_per_nm_;
  public:
    cdom(double abs_coef_440 = 0.01, double slope_per_nm = 0.017)
      : abs_coef_440_{abs_coef_440}, slope_per_nm_{slope_per_nm} {
    }
    double absorption_coefficient() const {
      double wl_nm = wavelength()*1e9;
      return abs_coef_440_ * exp(-slope_per_nm_*(wl_nm-440));
    }
    double scattering_coefficient() const {
      return 0;
    }
    double real_refractive_index() const {
      return 1.0;
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const
    {
      return rayleigh_mueller(angle(scattering_direction), 0.039);
    }
  };
}
}

#endif
