#ifndef flick_radiator_toa_solar
#define flick_radiator_toa_solar

#include "radiator.hpp"

namespace flick {
  namespace radiator {
    class toa_solar : public radiator {
      const std::string path_{"/radiator"};
    public:
      toa_solar() {
	spectrum_ = read<pp_function>(path_+"/toa_solar.txt");
      }
      double irradiance(double wavelength) {
	return spectrum_.value(wavelength);
      }
    };
  }
}

#endif

