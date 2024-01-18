#ifndef flick_delta_fit
#define flick_delta_fit

#include <numbers>
#include "../linalg/solve_with_eigen.hpp"
#include "../wigner/wigner_fit.hpp"
#include "legendre.hpp"
 
namespace flick {
  template<class Function>
  class delta_fit {
    std::vector<double> coefficients_;
  public:
    delta_fit(const Function& f, int n_terms)
      : coefficients_(n_terms) {
      int n_angles = wigner_n_sampling_points(n_terms);
      double forward_max = wigner_mu_max(n_angles);
      std::vector<double> x = range(-1,forward_max,n_angles).linspace();
      linalg::matrix m(x.size(), std::vector<double>(n_terms));
      legendre p(n_terms, x);
      for (size_t i=0; i<x.size(); ++i) {
	for (size_t j=0; j<n_terms; ++j) {
	  m[i][j] = p.value(j,i)/f.value(x[i]); 
	}
      }
      auto v = std::vector<double>(x.size(),1);
      coefficients_ = linalg::solve(m,v);
    }
    std::vector<double> coefficients() {
      return coefficients_;
    }
    std::vector<double> function_values(const std::vector<double>& x) {
      return legendre_evaluation(coefficients_).values(x);
    }
    double scaling_factor() const {
      return coefficients_.at(0)*4*std::numbers::pi;
    }
  };  
}

#endif
