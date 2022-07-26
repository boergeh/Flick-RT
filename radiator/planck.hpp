#ifndef flick_radiator_planck
#define flick_radiator_planck

#include "../numeric/constants.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"
#include "../environment/exception.hpp"

namespace flick {
   class planck
  // Black-body irradiance spectrum
  {
    double T_{5800};
    pl_function inv_acc_;
  public:
    planck(double temperature)
      : T_{temperature} {}
    double irradiance(double wavelength) {
      double wl = wavelength;
      using namespace constants;
      return 2*pi*h*pow(c,2)/pow(wl,5)/(exp(h*c/(wl*k_B*T_))-1);
    }
    double inverted_accumulated_wavelength(double fraction) {
      const size_t n_points{300};
      if (inv_acc_.size() < n_points) {
	auto x = range(0.01,20,n_points).linspace();
	auto y = std::vector<double>(n_points);
	for (size_t i=0; i < n_points; ++i) {
	  y[i] = pow(x[i],3)/(exp(x[i])-1);
	}
	auto f = pl_function{x,y};
	inv_acc_ = inverted_accumulation(f);
      }
      double x = inv_acc_.value((1-fraction)*total_integral());
      using namespace constants;
      return h*c/(x*k_B*T_);
    }
    pp_function irradiance_spectrum(size_t n_points) {
      if (n_points < 2)
	throw exception("radiator_planck");
      std::vector<double> wls(n_points);
      std::vector<double> irr(n_points);
      std::vector<double> area_fraction = range(0,1,n_points+1).linspace();
      double da = (area_fraction[1]-area_fraction[0])/2;
      for (size_t i=0; i < n_points; ++i) {
	wls[i] = inverted_accumulated_wavelength(area_fraction[i]+da);
	irr[i] = irradiance(wls[i]);
      }
      return pp_function{wls,irr};
    }
  private:
    double total_integral()
    // Using https://en.wikipedia.org/wiki/Riemann_zeta_function with
    // s=4 to integrage x^3/(e^x-1).
    {
      return pow(constants::pi,4)/15;
    }
  };
}

#endif

