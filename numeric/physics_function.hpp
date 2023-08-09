#ifndef flick_physics_functions
#define flick_physics_functions

#include "constants.hpp"

namespace flick {
  class henyey_greenstein
  // Henyey-Greenstein scattering phase function
  {
    double asymmetry_factor_{0};
  public:
    henyey_greenstein(double g) : asymmetry_factor_{g} {}
    double phase_function(double theta) const
    {
      return value(cos(theta)); 
    }
    double value(double mu) const
    // Note that integral over 4*pi equals one 
    {
      double g = asymmetry_factor_;
      double arg = 1+pow(g,2)-2*g*mu;
      return 1/(4*constants::pi)*(1-pow(g,2))/pow(arg,3./2);
    }
    double inverted_accumulated_angle(double fraction) const 
    // Note that fraction has range [0 1], and range value zero
    // corresponds to angle zero
    {
      double g = asymmetry_factor_;
      if (fabs(g) < 1e-9) // isotropic scattering
	return acos(1-2*fraction);
      double arg = (1-pow(g,2))/(1-g+2*g*(1-fraction));
      double mu = (1+pow(g,2)-pow(arg,2))/(2*g);
      return acos(mu); 
    }
  };
  
  class size_distribution
  // Size number distribution. Integral over all sizes equals one.
  {
  protected:
    const double pi_{constants::pi};
    double a_;
    double b_;
  public:
    size_distribution(double a, double b)
      : a_{a}, b_{b} {} 
    virtual double center() const = 0;
    virtual double width() const = 0;
    virtual double value(double x) const = 0;
    virtual double weighted_integral(double alpha) const = 0;
    double particles_per_volume(double volume_fraction) const {
      return volume_fraction * 1/(4./3*pi_*weighted_integral(3));
    }
    std::vector<double> values(const std::vector<double>& xv) {
      std::vector<double> yv(xv.size());
      for (size_t i=0; i<yv.size(); i++) {
	yv[i] = value(xv[i]);
      }
      return yv;
    }
  };
  
  class log_normal_distribution : public size_distribution
  // https://en.wikipedia.org/wiki/Log-normal_distribution
  {
  public:
    log_normal_distribution(double mu, double sigma)
      : size_distribution(mu,sigma) {
    }
    double value(double x) const {
      return 1/(x*b_*sqrt(2*pi_))*exp(-pow(log(x)-a_,2)/(2*pow(b_,2))); 
    }
    double center() const {
      return exp(a_+3*pow(b_,2)); // median of size volume distribution
    }
    double width() const {
      return b_;
    }
    double weighted_integral(double alpha) const {
      return exp(a_*alpha+0.5*pow(alpha*b_,2));
    }
    double average_area() {
      return constants::pi*weighted_integral(2);
    }
    static std::tuple<double, double>
    from_volume_distribution(double mu, double sigma) {
      return {mu-3*pow(sigma,2), sigma};
    }
  };
}

#endif
