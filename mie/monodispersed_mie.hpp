#ifndef flick_monodispersed_mie
#define flick_monodispersed_mie

#include "mie.hpp"

namespace flick
  // Implementation based on the following two papers: (1) Mishchenko,
  // M.I. and Yang, P., 2018. Far-field Lorenz–Mie scattering in an
  // absorbing host medium: theoretical formalism and FORTRAN
  // program. Journal of Quantitative Spectroscopy and Radiative
  // Transfer, 205, pp.241-252. (2) Mishchenko, M.I., Dlugach, J.M.,
  // Lock, J.A. and Yurkin, M.A., 2018. Far-field Lorenz–Mie
  // scattering in an absorbing host medium. II: Improved stability of
  // the numerical algorithm. Journal of Quantitative Spectroscopy and
  // Radiative Transfer, 217, pp.274-277.
{
  class special_function {
  protected:
    stdcomplex z_;
    stdvectorc f_;
  public:
    special_function(const stdcomplex& z, int n_terms)
      : z_{z}, f_(n_terms) {
      assert(f_.size()>1);
    }
    stdvectorc terms() {
      return f_;
    }
    stdvectorc times_z_derivatives() {
      stdvectorc d(f_.size());
      d[0] = 0;
      for (size_t n=1; n<f_.size(); ++n) {
	d[n] = z_ * f_[n-1] - double(n) * f_[n];
      }
      return d;
    }
  };

  class spherical_hankel : public special_function {
  public:
    spherical_hankel(const stdcomplex& z, int n_terms)
      : special_function(z,n_terms) {
      f_[0] = -1i*exp(1i*z)/z;
      f_[1] = -exp(1i*z)*(z+1i)/pow(z,2);
      for (size_t n=1; n<f_.size()-1; ++n) {
	f_[n+1] = (2*n+1.)/z*f_[n]-f_[n-1];
      }
    }
  };

  class bessel_first_kind : public special_function {    
    stdcomplex r_asymptotic(const stdcomplex& z, size_t n_terms) {
      return z/(2*n_terms+1.);
    }
    stdvectorc r(const stdcomplex& z, int n_high, int n_low,
		 const stdcomplex& r_high) {
      stdvectorc v(n_high-n_low);
      v.end()[-1] = r_high;
      for (int n=v.size()-2; n>=0; --n) {
	v[n] = 1./((2*n+1.)/z - v[n+1]);
      }
      return v;
    }
  public:
    bessel_first_kind(const stdcomplex& z, int n_terms)
      : special_function(z,n_terms) {
      double error = 1;
      int n_extra = 2;
      stdcomplex r_previous = r_asymptotic(z,n_terms);
      while(log10(error) > -std::numeric_limits<double>::digits10) {
	stdcomplex r0 = r_asymptotic(z,n_terms+n_extra);
	stdcomplex r_extra = r(z,n_terms+n_extra,n_terms,r0).front();
	error = std::abs(r_extra/r_previous-1.0);
	n_extra *= 2;
	r_previous = r_extra;
      }
      stdvectorc r_stable = r(z,n_terms,0,r_previous);
      //std::cout << "r_prev: " << r_previous <<std::endl;
      //std::cout << "r_stable: " << r_stable <<std::endl;
      //std::cout << "dn = "<< n_extra << ", error = " << error << std::endl;

      f_[0] = sin(z)/z;
      for (size_t n=1; n<f_.size(); ++n) {
	f_[n] = r_stable[n] * f_[n-1];
      }
    }
  };

  class monodispersed_mie : public basic_monodispersed_mie {
    
    int n_terms_ = n_terms();

    int n_terms() {
      stdcomplex x = size_parameter_in_host();
      return std::abs(x) + 4.05*std::pow(std::abs(x),1./3) + 14.;
    }
    std::tuple<stdvector,stdvector> pi_tau_polynomials(double angle) {
      stdvector pi(n_terms_);
      stdvector tau(n_terms_+1);
      double u = cos(angle);
      pi[0] = 1;
      pi[1] = 3 * u;
      tau[0] = u;
      for (size_t n=1; n<n_terms_; ++n) {
	double s = u * pi[n];
	double t = s - pi[n-1];
	pi[n+1] = s + (n+1) / n * t;
	tau[n] = n * t - pi[n-1];
      }
      tau.pop_back();
      return {pi, tau};      
    }
  public:
    //using basic_monodispersed_mie::basic_monodispersed_mie;
    
    monodispersed_mie(const stdcomplex& m_host,
		      const stdcomplex& m_sphere,
		      double vacuum_wl) :
      basic_monodispersed_mie::basic_monodispersed_mie(m_host,
						       m_sphere,
						       vacuum_wl) {
    }
    
    std::tuple<stdvectorc,stdvectorc> ab_coefficients() {
      //std::cout << "n_terms = " << n_terms_ << std::endl;
      //std::cout << "x: " << size_parameter_in_host() <<std::endl;

      stdcomplex m = m_sphere_ / m_host_;
      stdcomplex x = size_parameter_in_host();
      
      bessel_first_kind jx(x,n_terms_);
      stdvectorc jx_t = jx.terms();
      stdvectorc jx_d = jx.times_z_derivatives();
      
      bessel_first_kind jmx(m*x,n_terms_);  
      stdvectorc jmx_t = jmx.terms();
      stdvectorc jmx_d = jmx.times_z_derivatives();
      
      spherical_hankel hx(x,n_terms_);
      stdvectorc hx_t = hx.terms();
      stdvectorc hx_d = hx.times_z_derivatives();
      stdvectorc A = jmx_t * jx_d;
      stdvectorc B = jx_t * jmx_d;
      stdvectorc C = jmx_t * hx_d;
      stdvectorc D = hx_t * jmx_d;

      return {(pow(m,2)*A-B)/(pow(m,2)*C-D), (A-B)/(C-D)};      
    }
    void radius(double r) {
      radius_ = r;
      n_terms_ = n_terms();
    }
    void angles(const stdvector& angles) {
      angles_ = angles;
    }
    double absorption_cross_section() const {
      return 0;
    }
    double scattering_cross_section() const {
      return 0;
    }
    stdvector scattering_matrix_element(size_t row, size_t col) const
    // Note that integratinig element 0,0 over all 4*pi solid angles gives
    // the scattering cross section, where we count from 0 instead of one.
    {
      return {0};
    }    
  };
}
/*
function [a, b] = expansion_coefficients(m,x1,n_max)
  for n=1:n_max
    jmx = besselj(n+0.5,m*x1);
    jx = besselj(n+0.5,x1);
    hx = besselh(n+0.5,1,x1);
    dmxjmx = m*x1*besselj(n-1+0.5,m*x1)-n*besselj(n+0.5,m*x1);
    dxjx = x1*besselj(n-1+0.5,x1)-n*besselj(n+0.5,x1);
    dxhx = x1*besselh(n-1+0.5,1,x1)-n*besselh(n+0.5,1,x1);  
    A = jmx*dxjx;
    B = jx*dmxjmx;
    C = jmx*dxhx;
    D = hx*dmxjmx;
    a(n) = (m^2*A-B)/(m^2*C-D);
    b(n) = (A-B)/(C-D);
  end
end

function [pi, tau] = pi_tau_polynomials(theta, n_max)
  pi(1) = 1;
  pi(2) = 3*cos(theta);
  tau(1) = cos(theta);
  for n=2:n_max
    s = cos(theta)*pi(n);
    t = s-pi(n-1);
    pi(n+1) = s + (n+1)/n*t;
    tau(n) = n*t-pi(n-1);
  end
  pi(end)=[];
end

*/



#endif


