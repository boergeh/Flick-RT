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
}

#endif
