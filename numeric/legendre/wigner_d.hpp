#ifndef flick_wigner_d
#define flick_wigner_d

#include <vector>
//#include "../std_operators.hpp"

namespace flick {
  class wigner_d
  /* Weigner-d function, see appendix F.4 in: Mishchenko, M.I.,
     Travis, L.D. and Lacis, A.A., 2006. Multiple scattering of light
     by particles: radiative transfer and coherent
     backscattering. Cambridge University Press.
     https://www.giss.nasa.gov/staff/mmishchenko/publications/Book3.pdf
   */
  {
    stdvector t_;
    double x_;
    int m_, n_;
  public:
    wigner_d(double x, int m, int n, int n_terms)
      : x_{x}, m_{m}, n_{n}, t_(n_terms,0) {     
      int s_min = wigner_d::leading_zeros(m_,n_); 
      t_[s_min] = at_start_term(s_min);
      for (size_t s=s_min; s < t_.size()-1; ++s) {
	double r1 = sqrt(pow(s+1,2)-pow(m_,2));
	double r2 = sqrt(pow(s+1,2)-pow(n_,2));
	double r3 = sqrt(pow(s,2)-pow(m_,2));
	double r4 = sqrt(pow(s,2)-pow(n_,2));
	double k1 = 1/(s*r1*r2);
	double k2 = (2*s+1)*(s*(s+1)*x_-m_*n_);
	double k3 = (s+1)*r3*r4;
	if (s==0)
	  t_[s+1] = x_;
        else
	  t_[s+1] = k1*(k2*t_[s] - k3*t_[s-1]);
      }
    }
    stdvector terms() {
      return t_;
    }
    static size_t leading_zeros(int m, int n) {
      return std::max(abs(m),abs(n));
    }
  private:
    double at_start_term(int s) {
      int xi = 1;
      if (n_ < m_)
	xi = pow(-1,m_-n_);	
      int a = abs(m_-n_);
      int b = abs(m_+n_);
      return xi * pow(2,-s)
	* pow(factorial(2*s)/factorial(a)/factorial(b),0.5)
	* pow(1-x_,0.5*a) * pow(1+x_,0.5*b);
    }
    double factorial(int n) {
      return std::tgamma(n+1);
    }
  };
}

#endif
