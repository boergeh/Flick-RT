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
      if (mu >1 or mu<-1)
	throw std::out_of_range("henyey greenstein function mu out of range");
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

  class fournier_forand
  // Fournier, G.R. and Forand, J.L., 1994, October. Analytic phase
  // function for ocean water. In Ocean Optics XII (Vol. 2258,
  // pp. 194-201). SPIE.
  // See also Ocean Optics Web Book
  {
    double junge_slope_;
    double refractive_index_; // real refractive index relative to water
  public:
    fournier_forand(double junge_slope, double refractive_index)
      : junge_slope_{junge_slope}, refractive_index_{refractive_index} {
    }
    fournier_forand(double asymmetry_parameter)
    // Fournier, G.R., 2011. Derivation of an explicit expression for
    // the Fournier-Forand phase function in terms of the mean cosine.
    {
      double g = asymmetry_parameter;
      ensure((1-g) > 0.0002 and (1-g) < 0.6);
      double a = 1-g;
      double b = a/(23-7.5*a);
      refractive_index_ = pow(b,2./5) + 1;
      junge_slope_ = 6*(refractive_index_-1)+3;
    }
    double value(double mu) const {
      using namespace constants;      
      double d = delta(mu);
      double n = nu();
      double a = 1-d;
      double b = pow(d,n);
      double c = 1-b;
      double d180 = delta(-1);
      double b180 = pow(d180,n);
      double k1 = 1/(4*pi*pow(a,2)*b);
      double k2 = (1-b180)/(16*pi*(d180-1)*b180);
      double p = k1*(n*a-c+(d*c-n*a)/sin_squared_factor(mu)) + k2*(3*pow(mu,2)-1);
      return p;
    }
  private:
    void ensure(bool b) const {
      if (!b)
	throw std::runtime_error("fournier_forand");
    }
    double nu() const {
      return  (3 - junge_slope_)/2;
    }
    double delta(double mu) const {
      return 4 / (3*pow(refractive_index_-1,2)) * sin_squared_factor(mu);
    }
    double sin_squared_factor(double mu) const {
      return (1-mu)/2;
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
