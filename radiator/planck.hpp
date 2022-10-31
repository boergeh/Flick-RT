#ifndef flick_radiator_planck
#define flick_radiator_planck

#include "radiator.hpp"

namespace flick {
  namespace radiator {
    class planck : public radiator
    // Black-body irradiance spectrum
    {
      double T_{5800};
    public:
      planck(double temperature) : T_{temperature} {
	size_t n_points = 500;
	auto x = range(0.01,20,n_points).linspace();
	auto y = std::vector<double>(n_points);
	for (size_t i=0; i < n_points; ++i) {
	  y[i] = pow(x[i],5)/(exp(x[i])-1);
	}
	auto wl = std::vector<double>(n_points);
	auto pl = std::vector<double>(n_points);
	using namespace constants;
	double kx = h*c/(k_B*T_);
	double ky = 2*pi*pow(k_B*T_,5)/(pow(h,4)*pow(c,3));
	for (size_t i=0; i < n_points; ++i) {
	  wl[i] = kx / x[n_points-i-1];
	  pl[i] = ky * y[n_points-i-1];
	}
	spectrum_ = pp_function{wl,pl};
      }
    };
  }
}

#endif

