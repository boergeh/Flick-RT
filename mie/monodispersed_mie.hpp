#ifndef flick_monodispersed_mie
#define flick_monodispersed_mie

#include "basic_monodispersed_mie.hpp"
#include "../numeric/legendre/legendre.hpp"

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
      if (f_.size()<=1)
	throw std::runtime_error("monodispersed_mie");
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

  class spherical_bessel : public special_function {    
    stdcomplex r_asymptotic(const stdcomplex& z, size_t large_n) {
      return z/(2*large_n+1.);
    }
        
    stdvectorc r(const stdcomplex& z, size_t n_terms_high, int n_terms_low) {
      stdvectorc v(n_terms_high);
      v.end()[-1] = r_asymptotic(z,n_terms_high-1);
      for (int n=n_terms_high-2; n>=n_terms_low; --n) {
	v[n] = 1./((2*n+1.)/z - v[n+1]);
      }
      return v;
    }
  public:
    spherical_bessel(const stdcomplex& z, int n_terms)
      : special_function(z,n_terms) {
      double error = 1;
      int n_extra = 2;
      stdcomplex r_previous = r_asymptotic(z,n_terms);
      while(log10(error) > -(std::numeric_limits<double>::digits10)) {
	stdcomplex r_high = r(z,n_terms+n_extra,n_terms-1)[n_terms-1];
	error = std::abs(r_high/r_previous-1.0);
	n_extra *= 2;
	r_previous = r_high;
      }
      stdvectorc r_stable = r(z,n_terms+n_extra,2);
      f_[0] = sin(z)/z;
      f_[1] = sin(z)/pow(z,2)-cos(z)/z; // Avoids instability when sin(z)=0
      for (size_t n=2; n<f_.size(); ++n) {
	f_[n] = r_stable[n] * f_[n-1];
      }
    }
  };

  class monodispersed_mie : public basic_monodispersed_mie {  
    const double pi_ = constants::pi;

    int n_terms_;
    stdvectorc a_;
    stdvectorc b_;
    double C_ext_;
    double C_scat_;
    mutable stdvectorc S11_;
    mutable stdvectorc S22_;
    mutable bool changed_radius_ = true;
    mutable bool changed_angles_ = true;
    double refractive_index_slope_{0.0};
    
    std::tuple<stdvector,stdvector> pi_tau_polynomials(double angle) const {
      stdvector pi(n_terms_+1);
      stdvector tau(n_terms_);
      double u = std::cos(angle);
      pi[0] = 0; 
      pi[1] = 1;
      tau[0] = u;
      for (size_t n=1; n<n_terms_; ++n) {
	double s = u * pi[n];
	double t = s - pi[n-1];
	pi[n+1] = s + (n+1.) / n * t;
	tau[n] = n * t - pi[n-1];
      }
      pi.pop_back();
      return {pi, tau};      
    }
    void set_n_terms() {
      stdcomplex x = size_parameter_in_host();
      n_terms_ = 8. + std::abs(x) + 4.05*std::pow(std::abs(x),1./3);
    }
  public:
    using basic_monodispersed_mie::basic_monodispersed_mie;
 
    std::tuple<stdvectorc,stdvectorc> ab_coefficients() {
      stdcomplex m = m_sphere_ / m_host_;
      stdcomplex x = size_parameter_in_host();
      
      spherical_bessel jx(x,n_terms_);
      stdvectorc jx_t = jx.terms();
      stdvectorc jx_d = jx.times_z_derivatives();
      
      spherical_bessel jmx(m*x,n_terms_);  
      stdvectorc jmx_t = jmx.terms();
      stdvectorc jmx_d = jmx.times_z_derivatives();
      
      spherical_hankel hx(x,n_terms_);
      stdvectorc hx_t = hx.terms();
      stdvectorc hx_d = hx.times_z_derivatives();
      stdvectorc A = jmx_t * jx_d;
      stdvectorc B = jx_t * jmx_d;
      stdvectorc C = jmx_t * hx_d;
      stdvectorc D = hx_t * jmx_d;

      stdvectorc a = (pow(m,2)*A-B)/(pow(m,2)*C-D);
      stdvectorc b = (A-B)/(C-D);
      
      return {a, b};      
    }
    std::tuple<stdvectorc,stdvectorc> s_functions() const {
      stdvectorc S11(angles_.size(),stdcomplex{0,0});
      stdvectorc S22(angles_.size(),stdcomplex{0,0});
      for (size_t i=0; i<angles_.size(); ++i) {
	auto [pi, tau] = pi_tau_polynomials(angles_[i]);
	for (size_t n=1; n<n_terms_; ++n) {
	  stdcomplex c = 1i/wavenumber_in_host_*(2*n+1.)/(n*(n+1.));
	  S11[i] += c*(a_[n]*tau[n] + b_[n]*pi[n]);
	  S22[i] += c*(a_[n]*pi[n] + b_[n]*tau[n]);
	}
      }
      return {S11, S22};
    }
    std::tuple<double,double> es_coefficients() const {
      stdcomplex ext{0,0};
      double scat{0};
      const stdcomplex& k = wavenumber_in_host_;
      for (size_t n=1; n<n_terms_; ++n) {
	ext += (2*n+1.) * (a_[n] + b_[n]);
	scat += (2*n+1) * (norm(a_[n]) + norm(b_[n]));
      }
      double C_ext = 2*pi_/real(k)*real(ext/k); 
      double C_scat = 2*pi_/norm(k)*scat; 
      return {C_ext, C_scat};
    }
    void refractive_index_slope(double s) {
      refractive_index_slope_ = s;
    }
    void radius(double r) {
      radius_ = r;
      double r0_ = 1e-6;
      m_sphere_.real(m_sphere_at_r0_.real()*pow(radius_/r0_,refractive_index_slope_));
      set_n_terms();
      std::tie(a_, b_) = ab_coefficients();
      std::tie(C_ext_, C_scat_) = es_coefficients();
      changed_radius_ = true;
    }
    void angles(const stdvector& angles) {
      angles_ = angles;
      changed_angles_ = true;
    }
    void quadrature_angles(size_t n) {
      stdvector x = read_quadrature(n).column(0);
      std::reverse(x.begin(),x.end());
      angles(vec::acos(x));
    }
    stdvector angles() {
      return angles_;
    }
    double absorption_cross_section() const {
      return C_ext_ - C_scat_;
    }
    double scattering_cross_section() const {
      return C_scat_;
    }
    stdvector scattering_matrix_element(size_t row, size_t col) const
    // Note that integratinig element F11 over all 4*pi solid angles gives
    // the scattering cross section.
    {
      if (changed_angles_ or changed_radius_) {
	std::tie(S11_, S22_) = s_functions();
	changed_radius_ = false;
	changed_angles_ = false;
      }
      
      bool F11 = (row==0 and col==0);
      bool F22 = (row==1 and col==1);
      bool F33 = (row==2 and col==2);
      bool F44 = (row==3 and col==3);
      bool F12 = (row==0 and col==1);
      bool F21 = (row==1 and col==0);
      bool F34 = (row==2 and col==3);
      bool F43 = (row==3 and col==2);

      using namespace vec;
      if (F11 or F22)
	return 0.5*((vec::abs(S11_)^2)+(vec::abs(S22_)^2));
      else if (F33 or F44)
	return real(S11_*conj(S22_));
      else if (F12 or F21)
	return 0.5*((vec::abs(S11_)^2)-(vec::abs(S22_)^2));
      else if (F34)
	return imag(S11_*conj(S22_));
      else if (F43)
	return -1*imag(S11_*conj(S22_));
      else
	return stdvector(angles_.size(),0);
    }    
  };
}

#endif


