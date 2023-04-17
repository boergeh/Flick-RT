#ifndef flick_legendre
#define flick_legendre

#include <vector>
#include "../../environment/input_output.hpp"
#include "../function.hpp"
#include "../std_operators.hpp"
#include "wigner_d.hpp"

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
    two_columns quadrature_;
    Function& f_;
  public:
    gl_integral_vector(Function& f, size_t n_points)
      : f_{f} {
      quadrature_ = read_quadrature(n_points);
    }
    std::tuple<stdvector,std::vector<stdvector>>
    xy_integration_points(double from, double to) {
      stdvector x = from + (quadrature_.column(0)+1)/2*(to-from);
      std::vector<stdvector> y(f_.size(),stdvector(x.size()));
      for (size_t i = 0; i < x.size(); ++i) {
	stdvector f = f_.value(x[i]);
	for (size_t j = 0; j < f.size(); ++j) {
	   y[j][i] = f[j];
	}
      }
      return {x,y};
    }
    stdvector value(double from, double to) {
      auto [x,y] = xy_integration_points(from,to);
      stdvector v(y.size());
      for (size_t i = 0; i < v.size(); ++i) {
	v[i] = vec::sum(quadrature_.column(1)*y[i])/2*(to-from);
      }
      return v;   
    }
    stdvector of_abs_integrand(double from, double to) {
      auto [x,y] = xy_integration_points(from,to);
      stdvector v(y.size());
      for (size_t i = 0; i < v.size(); ++i) {
	v[i] = vec::sum(quadrature_.column(1)*vec::abs(y[i]))/2*(to-from);
      }
      return v;   
    }
  };
  
  class legendre {
    std::vector<stdvector> p_;
  public:
    legendre(size_t n_terms, stdvector x) {
      p_.resize(x.size());
      for (size_t i=0; i < x.size(); ++i) {
	p_[i] = wigner_d(x[i],0,0,n_terms).terms();
      }
    }
    double value(size_t term_number, size_t x_number) {
      return p_[x_number][term_number];
    }
  };

  template<class Function>
  stdvector legendre_expansion(const Function& f,
					 size_t n_terms,
					 size_t n_points = 64) {
    stdvector terms(n_terms);
    stdvector x = read_quadrature(n_points).column(0);
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
    const stdvector& coefficients_;
  public:
    legendre_evaluation(const stdvector& coefficients)
      : coefficients_{coefficients} {}
    stdvector values(const stdvector& x) {
      stdvector values(x.size());
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
}

#endif
