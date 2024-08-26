#ifndef flick_spherical_harmonics
#define flick_spherical_harmonics

#include "../wigner/wigner_d.hpp"
#include <complex>
#include <numbers>

namespace flick {
  struct direction {
    double theta;
    double phi;
  };
  class spherical_harmonics {
    std::vector<std::vector<std::complex<double>>> yml_;
    int l_max_;
  public:
    spherical_harmonics(const direction& d, int l_max)
      : l_max_{l_max} {
      std::complex<double> i(0.0,1.0);
      for (int m=-l_max ; m<=l_max; ++m) {
	std::vector<double> dml = wigner_d(cos(d.theta),m,0,l_max+1).terms();
	std::vector<std::complex<double>> dmlc(dml.size());
	for (size_t l=0; l<dml.size();++l) {
	  dmlc[l] = dml[l]*sqrt((2*l+1)/(4*std::numbers::pi))*exp(m*d.phi*i);
	}
	yml_.push_back(dmlc);
      }
    }
    std::complex<double> lm(size_t l, int m) {
      return yml_[m+l_max_][l];
    }
  };

  class associated_legendre {
    spherical_harmonics sh_;
  public:
    associated_legendre(double theta, int l_max)
      : sh_(direction{theta,0},l_max) {}
    double lm(size_t l, int m) {
      return std::real(sh_.lm(l,m)) / sqrt((2*l+1)*factorial(l-m)/(4*std::numbers::pi*factorial(l+m)));
    }
  private:
    double factorial(int n) {
      return std::tgamma(n+1);
    }
  };
  
  class associated_legendre_curve {
    std::vector<associated_legendre> al_;
  public:
    associated_legendre_curve(size_t n_theta, int l_max) {
      auto theta = linspace(0,std::numbers::pi,n_theta);
      for (size_t i=0; i<n_theta; ++i) {
	al_.push_back(associated_legendre(theta[i],l_max));
      }
    }
    std::vector<double> lm(size_t l, int m) {
      std::vector<double> v;
      for (size_t i=0; i<al_.size(); ++i) {
	v.push_back(al_[i].lm(l,m));
      }
      return v;
    }
  private:
    double factorial(int n) {
      return std::tgamma(n+1);
    }
    std::vector<double> linspace(double from, double to, size_t n) {
      double step = (to-from)/n;
      std::vector<double> v(n);
      v[0] = from;
      for (size_t i=1; i<v.size(); i++)
	v[i] = v[i-1] + step;
      return v;
    }
  };
}

#endif
