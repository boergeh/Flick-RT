#ifndef flick_sh_expansion
#define flick_sh_expansion

#include "spherical_harmonics.hpp"

namespace flick {
  class sh_expansion {
    size_t l_max_;
    using complex_vector = std::vector<std::complex<double>>;
    std::vector<complex_vector> coefficients_;
    double pi_ = std::numbers::pi;
  public:
    sh_expansion(size_t l_max)
      : l_max_{l_max}, coefficients_(l_max+1,complex_vector(2*l_max+1,std::complex<double>(0,0))) {
    }
    void set_coefficient_lm(size_t l, int m, const std::complex<double>& c) {
      coefficients_[l][l_max_+m] = c;
    }
    double value(const direction& d) {
      spherical_harmonics y(d,coefficients_.size()-1);
      double v = 0;
      for (size_t l=0; l<coefficients_.size(); ++l) {
	for (size_t n=0; n<coefficients_[l].size(); ++n) {
	  int m = n-l_max_;
	  v += std::abs(coefficients_[l][n] * y.lm(l,m));
	}
      }
      return v;
    }
    std::vector<double> theta_distribution(size_t n) {
      double d_theta = pi_/n/2;
      return linspace(d_theta, pi_-d_theta, n);
    }
    std::vector<double> phi_distribution(size_t n) {
      double d_phi = pi_/n;
      return linspace(d_phi, 2*pi_-d_phi, n);
    }
    std::vector<std::vector<double>> value_distribution(size_t n_theta, size_t n_phi) {
      std::vector<std::vector<double>> r(n_theta, std::vector<double>(n_phi));
      auto theta = theta_distribution(n_theta);
      auto phi = phi_distribution(n_phi);     
      for (size_t i=0; i<n_theta; ++i) {
	for (size_t j=0; j<n_phi; ++j) {
	  direction d{theta[i], phi[j]};
	  r[i][j] = value(direction{theta[i], phi[j]});
	}
      }
      return r;
    }
    void write_distribution(size_t n_theta, size_t n_phi, std::ostream& os, size_t precision=6) {
      os << std::setprecision(precision);
      auto theta = theta_distribution(n_theta);
      auto phi = phi_distribution(n_phi);
      auto v = value_distribution(n_theta, n_phi);
      os << 0 << " ";
      for (size_t j = 0; j < n_phi; ++j) {
	os <<  phi[j] << " ";
      }
      os << '\n';
      for (size_t i = 0; i < n_theta; ++i) {
	os << theta[i] << " ";
	for (size_t j = 0; j < n_phi; ++j) {
	  os << v[i][j] << " ";
	}
	os << '\n';
      }
    }
  private:
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
