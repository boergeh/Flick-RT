#ifndef flick_material_pure_water
#define flick_material_pure_water

#include "../material.hpp"
#include "../../polarization/rayleigh_mueller.hpp"

namespace flick {
namespace material {
  class pure_ice : public base {
    pp_function absorption_coefficient_;
    pp_function real_refractive_index_;
    const std::string path_{"/material/ice"};
  public:
    pure_ice() {
     absorption_coefficient_ = read<pp_function>
       (path_+"/absorption_coefficient.txt"); 
     real_refractive_index_ = read<pp_function>
       (path_+"/refractive_index.txt"); 
     absorption_coefficient_.add_extrapolation_points(1);
     real_refractive_index_.add_extrapolation_points(1);
    }
    double absorption_coefficient() const {
      return absorption_coefficient_.value(wavelength());
    }
    double scattering_coefficient() const {
      return pow(129.0 / (wavelength() * 1e9), 4.32);
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      return rayleigh_mueller(angle(scattering_direction), 0.039);
    }
    double real_refractive_index() const {
      return real_refractive_index_.value(wavelength());
    }
  };
}
}

#endif
