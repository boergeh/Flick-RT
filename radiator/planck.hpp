#ifndef flick_planck
#define flick_planck

#include "../numeric/constants.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"

namespace flick {
  class planck {
    double T_{6000};
    pl_function inv_acc_;
    size_t points_{200};
  public:
    planck(double temperature)
      : T_{temperature} {
      auto x = range(0,20,points_).linspace();
      auto y = std::vector<double>(points_);
      for (size_t i=0; i < points_; ++i) {
	y[i] = f(x[i]);
      }
      pl_function f{x,y};
      //inv_acc_ = inverted_accumulation(f); 

    }
    double value(double wavelength) {
      using namespace constants;
      double wl = wavelength;
      return 2*pi*h*pow(c,2)/pow(wl,5)/(exp(h*c/(wl*k_B*T_))-1);
    }
    double wavelength(double integral_fraction) {
      using namespace constants;
      double x = inv_acc_.value(integral_fraction);
      return h*c/(x*k_B*T_);
    }
    double f(double x) {
      return pow(x,3)/(exp(x)-1);
    }
    /*
    double x_max(double epsilon) {
      double df = 100;
      double x = 1;
      total_area = pow(pi,4)/15;
      while(df > epsilon) {
	df = std::rieman_zetha(4)*std::tgamma(4)/total_area;
	x += 1;
      }
      return x;	
    }
    */
  };
}

#endif
