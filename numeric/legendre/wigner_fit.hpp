#ifndef flick_wigner_fit
#define flick_wigner_fit

#include "../std_operators.hpp"
#include "wigner_d.hpp"
#include <armadillo>
 
namespace flick {
  enum class fit {
    absolute,
    relative,
    ratio_function
  };
  template<class Function>
  class wigner_fit {
    stdvector coefficients_;
    int m_, n_;
    pe_function ratio_function_;//{{-1,1},{1,1}};
  public:
    wigner_fit(const Function& f, int m, int n, int n_terms, fit fit,
	       const pe_function& rf = pe_function{{-1,1},{1,1}})
      : m_{m}, n_{n}, coefficients_(n_terms) {
      ratio_function_ = rf;
      size_t n_points = wigner_fit::n_sample_points(n_terms);
      stdvector x = range(-1,1,n_points).linspace();   
      arma::mat matrix(x.size(), n_terms);
      arma::vec v(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	stdvector d = wigner_d(x[i], m_, n_, n_terms).terms();
	if (fit == fit::absolute)
	  v(i) = f.value(x[i]);
	else if (fit == fit::relative)
	  v(i) = 1;
	else if (fit == fit::ratio_function)
	  v(i) = f.value(x[i])/ratio_function_.value(x[i]);
	  
	for (size_t j=0; j<n_terms; ++j) {
	  if (fit == fit::absolute)
	    matrix(i,j) = d[j];
	  else if (fit == fit::relative)
	    matrix(i,j) = d[j]/f.value(x[i]);
	  else if (fit == fit::ratio_function)
	    matrix(i,j) = d[j]/ratio_function_.value(x[i]);
	}
      }
      size_t n_zeros = wigner_d::leading_zeros(m_,n_);
      if (n_zeros > 0)
	matrix.shed_cols(0,n_zeros-1);
      arma::vec c = arma::solve(matrix,v);
      if (n_zeros > 0)
	c.insert_rows(0,arma::vec(n_zeros,arma::fill::zeros));
      coefficients_ = arma::conv_to<std::vector<double>>::from(c);
    }
    void ratio_function(const pe_function& rf) {
      ratio_function_ = rf;
    }
    stdvector coefficients() {
      return coefficients_;
    }
    double value(double x) {
      stdvector d = wigner_d(x, m_, n_, coefficients_.size()).terms();
      return vec::sum(d*coefficients_);
    }
    static size_t n_sample_points(size_t n_terms) {
      return pow(n_terms,1.6);
    }
  };

  stdvector wigner_evaluate(const stdvector& coefficients, const stdvector& x, int m, int n) {
    stdvector v(x.size());
    for (size_t i=0; i<x.size(); ++i)
      v[i] = vec::sum(wigner_d(x[i],m,n,coefficients.size()).terms() * coefficients);
    return v;
  }
}

#endif
