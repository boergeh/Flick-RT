#ifndef flick_radiator_cie_a
#define flick_radiator_cie_a

#include "radiator.hpp"

namespace flick {
  namespace radiator {
    class cie_a : public radiator {
      const std::string path_{"/radiator"};
    public:
      cie_a() {
	spectrum_ = read<pp_function>(path_+"/cie_a.txt");
	spectrum_.scale_x(1e-9);
      }
      double irradiance(double wavelength) const {
	return spectrum_.value(wavelength);
      }
    };
  }
}

#endif

