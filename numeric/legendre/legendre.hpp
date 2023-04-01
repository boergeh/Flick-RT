#ifndef flick_legendre
#define flick_legendre

#include <vector>
#include "../../environment/input_output.hpp"
#include "../function.hpp"
#include "../std_operators.hpp"

namespace flick {
  two_columns read_quadrature(size_t n_points) {
    std::string fname ="numeric/legendre/gl_quadrature/q_"+
      std::to_string(n_points)+".txt";
    return read<two_columns>(fname);     
  }
  
  template<class Function>
  class gl_integral {
  protected:
    size_t n_points_;
    two_columns quadrature_;
    const Function& f_;
  public:
    gl_integral(const Function& f, size_t n_points)
      : f_{f}, n_points_{n_points} {
      quadrature_ = read_quadrature(n_points);
    }
    double value(double from, double to) const {
      double range = to - from;
      double v = 0;
      for (size_t n = 0; n < n_points_; ++n) {
	double x = quadrature_.column(0)[n];
	double w = quadrature_.column(1)[n];
	double local_x = from + (x+1)/2*range;
	v += w * f_.value(local_x);
      }
      return v/2*range;
    }
  };
  
  template<class Function>
  class gl_integral_vector {
  protected:
    size_t n_points_;
    two_columns quadrature_;
    Function& f_;
  public:
    gl_integral_vector(Function& f, size_t n_points)
      : f_{f}, n_points_{n_points} {
      quadrature_ = read_quadrature(n_points);
    }
    std::vector<double> value(double from, double to) {
      double range = to-from;
      const std::vector<double>& x0 = quadrature_.column(0);
      std::vector<double> x_shifted = from + (x0+1)/2*range;
      const std::vector<double>& weights = quadrature_.column(1);

      size_t n_out = f_.size();;
      std::vector<std::vector<double>> m;
      m.resize(n_out, std::vector<double>(n_points_));
      for (size_t i = 0; i < n_points_; ++i) {
	std::vector<double> y = f_.value(x_shifted[i]);
	for (size_t j = 0; j < y.size(); ++j) {
	   m[j][i] = weights[i] * y[j];
	}
      }
      std::vector<double> v(n_out);
      for (size_t j = 0; j < v.size(); ++j) {
	v[j] = vec::sum(m[j])/2*range;
      }
      return v;   
    }
  };
  
  class legendre {
    std::vector<std::vector<double>> p_;
  public:
    legendre(size_t n_terms, std::vector<double> x) {
      p_.resize(x.size());
      for (size_t i=0; i < x.size(); ++i) {
	std::vector<double> P(n_terms);
	P[0] = 1;
	P[1] = x[i];
	for (size_t l=1; l < n_terms-1; ++l)
	  P[l+1] = ((2.*l+1)*x[i]*P[l]-l*P[l-1])/(l+1.);
	p_[i] = P;
      }
    }
    double value(size_t term_number, size_t x_number) {
      return p_[x_number][term_number];
    }
  };

  template<class Function>
  std::vector<double> legendre_expansion(const Function& f,
					 size_t n_terms,
					 size_t n_points = 64) {
    std::vector<double> terms(n_terms);
    std::vector<double> x = read_quadrature(n_points).column(0);
    legendre legendre(n_terms, x);
    for (size_t i=0; i < terms.size(); ++i) {
      pl_function plf;
      for (size_t j=0; j < x.size(); ++j) {
	double y = f.value(x[j]) * legendre.value(i,j);
	plf.append(point{x[j],y});
      }
      terms[i] = (2.*i+1)/2 * gl_integral(plf,x.size()).value(-1,1);
    }
    return terms;
  }

  class legendre_evaluation {
    const std::vector<double>& coefficients_;
  public:
    legendre_evaluation(const std::vector<double>& coefficients)
      : coefficients_{coefficients} {}
    std::vector<double> values(const std::vector<double>& x) {
      std::vector<double> values(x.size());
      legendre legendre(coefficients_.size(), x);
      for (size_t i=0; i < x.size(); ++i) {
	double v = 0;
	for (size_t l=0; l < coefficients_.size(); ++l) {
	  v += coefficients_[l]*legendre.value(l,i);
	}
	values[i] = v;
      }
      return values;
    }
  };

  /*
  template<class Function>
  std::vector<double> delta_log(const Function& f,
				size_t n_terms,
				size_t n_points = 64) {
    std::vector<double> x = read_quadrature(n_points).column(0);
    pl_function pef;
    for (size_t i=0; i<x.size(); ++i) {
      pef.append({x[i],log(f.value(x[i]))});
    }
    std::vector<double> terms = legendre_expansion(pef,n_terms,n_points);
    std::vector<double> y = legendre_evaluation(terms).values(x);
    pef.clear();
    for (size_t i=0; i<x.size(); ++i) {
      pef.append({x[i],exp(y[i])});
    }
    return legendre_expansion(pef,n_terms,n_points);
  }
  */
}

#endif
