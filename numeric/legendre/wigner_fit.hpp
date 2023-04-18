#ifndef flick_wigner_fit
#define flick_wigner_fit

#include "../constants.hpp"
#include "wigner_d.hpp"
#include <armadillo>
 
namespace flick {
  enum class fit {
    absolute,
    relative
  };
  template<class Function>
  class wigner_fit {
    stdvector coefficients_;
    int m_, n_;
  public:
    wigner_fit(const Function& f, int m, int n, int n_terms, fit fit)
      : m_{m}, n_{n}, coefficients_(n_terms) {
      int n_points = pow(n_terms,1.6);
      stdvector x = range(-1,1,n_points).linspace();      
      arma::mat matrix(x.size(), n_terms);
      arma::vec v = arma::ones(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	stdvector d = wigner_d(x[i], m_, n_, n_terms).terms();
	if (fit == fit::absolute)
	  v(i) = f.value(x[i]);
	for (size_t j=0; j<n_terms; ++j) {
	  if (fit == fit::relative)
	    matrix(i,j) = d[j]/f.value(x[i]);
	  else
	    matrix(i,j) = d[j];
	}
      }       
      arma::vec c = arma::solve(matrix,v);
      coefficients_ = arma::conv_to<std::vector<double>>::from(c);
    }
    stdvector coefficients() {
      return coefficients_;
    }
    double value(double x) {
      stdvector d = wigner_d(x, m_, n_, coefficients_.size()).terms();
      return vec::sum(d*coefficients_);
    }
  };
}

#endif
