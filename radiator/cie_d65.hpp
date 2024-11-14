#ifndef flick_radiator_cie_d65
#define flick_radiator_cie_d65

#include "radiator.hpp"

namespace flick {
  namespace radiator {
    class cie_d65 : public radiator {
      const std::string path_{"/radiator"};
    public:
      cie_d65() {
	spectrum_ = read<pp_function>(path_+"/cie_d65.txt");
	spectrum_.scale_x(1e-9);
      }
      double irradiance(double wavelength) const {
	return spectrum_.value(wavelength);
      }
    };
  }
}

#endif

