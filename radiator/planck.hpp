#ifndef flick_planck
#define flick_planck

#include "../numeric/constants.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"
#include "../environment/exception.hpp"

namespace flick {
  // should make radiator base class
  
  class planck
  // Black-body irradiance spectrum
  {
    double T_{5800};
    pl_function inv_acc_;
    size_t n_points_{300};
    double x_min{0.01};
    double x_max{20};
  public:
    planck(double temperature)
      : T_{temperature} {
      auto x = range(x_min,x_max,n_points_).linspace();
      auto y = std::vector<double>(n_points_);
      for (size_t i=0; i < n_points_; ++i) {
	// move this into a general distribution in numeric
	y[i] = pow(x[i],3)/(exp(x[i])-1);
      }
      auto f = pl_function{x,y};
      inv_acc_ = inverted_accumulation(f);
    }
    double irradiance(double wavelength) {
      double wl = wavelength;
      using namespace constants;
      return 2*pi*h*pow(c,2)/pow(wl,5)/(exp(h*c/(wl*k_B*T_))-1);
    }
    double density_wavelength(double integral_fraction) {
      double x = inv_acc_.value((1-integral_fraction)*total_integral());
      using namespace constants;
      return h*c/(x*k_B*T_);
    }
    pp_function irradiance_spectrum(size_t n_points) {
      if (n_points < 2)
      	throw exception("planck radiator"); 
      std::vector<double> wls(n_points);
      std::vector<double> irr(n_points);
      std::vector<double> area_fraction = range(0,1,n_points+1).linspace();
      double da = (area_fraction[1]-area_fraction[0])/2;
      for (size_t i=0; i < n_points; ++i) {
	wls[i] = density_wavelength(area_fraction[i]+da);
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
