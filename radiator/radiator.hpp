#ifndef flick_radiator
#define flick_radiator

#include "../environment/input_output.hpp"
#include "../environment/exception.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"
#include "../numeric/constants.hpp"

namespace flick {
  namespace radiator {
    class radiator {
    protected:
      pp_function spectrum_;     
    public:
      double irradiance(double wavelength) const {
	return spectrum_.value(wavelength);
      }
      const pp_function& spectrum() const {
	return spectrum_;
      }
      pp_function spectrum(size_t n_points) const {	
	return importance_sampled(spectrum_, n_points);
      }
    };
  }
}

#endif

