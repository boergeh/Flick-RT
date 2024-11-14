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
      planck(double temperature)
	: T_{temperature} {
	initialize(0.01, 20, 500);
      }
      planck(double temperature, double wl_low, double wl_high, size_t n_points)
	: T_{temperature} {
	initialize(x(wl_high), x(wl_low), n_points);
      }
      double value(double wavelength) {
	return f(x(wavelength));
      }
    private:
      void initialize(double x_low, double x_high, size_t n_points) {
	std::vector<double> xv = range(x_low, x_high, n_points).linspace();
	for (int i=n_points-1; i >= 0; i--) {
	  point p(wavelength(xv[i]),f(xv[i]));
	  spectrum_.append(p);
	}
      }
      double kernel(double x) {
	return pow(x,5)/(exp(x)-1);
      }
      double k() {
	using namespace constants;
	return h*c/(k_B*T_);
      }
      double wavelength(double x) {
	return k()/x;
      }
      double x(double wavelength) {
	return k()/wavelength;
      }
      double f(double x) {
	using namespace constants;
	double kf = 2*pi*pow(k_B*T_,5)/(pow(h,4)*pow(c,3));
        return kf*kernel(x);
	
      }
    };
  }
}

#endif

