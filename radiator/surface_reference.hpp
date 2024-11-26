#ifndef flick_radiator_surface_reference
#define flick_radiator_surface_reference

#include "radiator.hpp"

namespace flick {
  namespace radiator {
    class surface_reference : public radiator {
      const std::string path_{"/radiator"};
    public:
      surface_reference() {
	spectrum_ = read<pp_function>(path_+"/surface_reference.txt");
	spectrum_.scale_x(1e-9);
	spectrum_.scale_y(1e9);
      }
      double irradiance(double wavelength) const {
	return spectrum_.value(wavelength);
      }
    };
  }
}

#endif

