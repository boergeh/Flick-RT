#ifndef flick_fresnel
#define flick_fresnel

#include "../numeric/vector.hpp"
#include "mueller.hpp"
#include <complex>
#include <cassert>

// See Wikipedia Fresnel equations. See also the fresnel_curves
// function below for further documentation.  Relatative refractive
// index is m = n_2/n_1, where radiation is incident from medium with
// refractive index n_1, and surface normal is pointing from medium
// with n_1 into medium with n_2.
namespace flick {
  using complex = std::complex<double>;
  class fresnel {
    std::complex<double> m_{1,0};
    double cos_theta_i_{1};
    complex cos_theta_t_{1,0};
  public:
    fresnel()=default;
    fresnel(const complex& relative_refractive_index, double cos_incidence_angle)
      : m_{relative_refractive_index}, cos_theta_i_{cos_incidence_angle}
    {
      if (cos_theta_i_ < -1)
	cos_theta_i_ = -1;
      else if (cos_theta_i_ > 1)
	cos_theta_i_ = 1;
      double x = 1-pow(cos_theta_i_,2);
      double sin_theta_i = sqrt(x);
      assert(sin_theta_i >= -1 && sin_theta_i <= 1);
      complex c = complex{sin_theta_i,0}/m_;
      cos_theta_t_ = sqrt(1.-c*c);
    }
    complex r_s() const {
      double a = cos_theta_i_;
      complex b = m_ * cos_theta_t_;
      return (a-b)/(a+b);
    }
    complex r_p() const {
      complex a = m_ * cos_theta_i_;
      const complex& b = cos_theta_t_;
      return (a-b)/(a+b);
    }
    complex t_s() const {
      return (2*cos_theta_i_)/(cos_theta_i_ + m_*cos_theta_t_);
    }
    complex t_p() const {
      return (2*cos_theta_i_)/(m_*cos_theta_i_ + cos_theta_t_);
    }
    double R_s() const {
      return std::pow(std::abs(r_s()),2);
    }
    double R_p() const {
      return std::pow(std::abs(r_p()),2);
    }   
    double R() const {
      return 0.5*(R_s() + R_p());
    }
    double T_s() const {
      return 1-R_s();
    }
    double T_p() const {
      return 1-R_p();
    }
    double T() const {
      return 0.5*(T_s() + T_p());
    }
  
    // Mueller matrix elements
    double R_11() const {return 0.5*(R_p()+R_s());}
    double R_12() const {return 0.5*(R_p()-R_s());}
    double R_33() const {return std::real(r_p()*std::conj(r_s()));}
    double R_34() const {return std::imag(r_p()*std::conj(r_s()));}
    double T_11() const {return 0.5*(T_p()+T_s());}
    double T_12() const {return 0.5*(T_p()-T_s());}
    double T_33() const {return std::real(t_p()*std::conj(t_s()));}
    double T_34() const {return std::imag(t_p()*std::conj(t_s()));} // check K-factor    
    double reflection_angle() const {
      return acos(cos_theta_i_);
    }
    complex transmission_angle() const {
      //assert(1-cos_theta_t_ > 0);
      return acos(cos_theta_t_);
    }
  };

  mueller reflection_mueller(const fresnel& f) {
    mueller m;
    m.add(0,0,f.R_11());
    m.add(0,1,f.R_12());
    m.add(1,0,f.R_12());
    m.add(1,1,f.R_11());
    m.add(2,2,f.R_33());
    m.add(2,3,f.R_34());
    m.add(3,2,-f.R_34());
    m.add(3,3,f.R_33());
    return m;
  }
  mueller transmission_mueller(const fresnel& f) {
    mueller m;
    m.add(0,0,f.T_11());
    m.add(0,1,f.T_12());
    m.add(1,0,f.T_12());
    m.add(1,1,f.T_11());
    m.add(2,2,f.T_33());
    m.add(2,3,f.T_34());
    m.add(3,2,-f.T_34());
    m.add(3,3,f.T_33());
    return m;
  }
  
  std::vector<double> linspace(double from, double to, size_t n) {
    std::vector<double> v(n);
    double step = (to-from)/(n-1);
    v[0] = from;
    for (size_t i = 1; i < n; ++i)
      v[i] = v[i-1] + step;
    return v;  
  }
  void fresnel_curves(std::complex<double> relative_refractive_index,
		      size_t n_angles, std::ostream& ostr) {
     using namespace constants;
    std::vector<double> theta = linspace(0,pi/2,n_angles);
    ostr << "# column 1: incidence angle in degrees" << std::endl
	 << "# column 2: s-component reflection amplitude coefficient (r_s)" << std::endl
	 << "# column 3: p-component reflection amplitude coefficient (r_p)" << std::endl
	 << "# column 4: s-component transmission amplitude coefficient (t_s)" << std::endl
	 << "# column 5: p-component transmission amplitude coefficient (t_p)" << std::endl
	 << "# column 6: s-component reflection power coefficient (R_s)" << std::endl
	 << "# column 7: p-component reflection power coefficient (R_p)" << std::endl
	 << "# column 8: s-component transmission power coefficient (T_s)" << std::endl
	 << "# column 9: p-component transmission power coefficient (T_p)" << std::endl;
    for (size_t i = 0; i < n_angles; ++i) {
      fresnel f(relative_refractive_index, theta.at(i));
      ostr << theta.at(i)/pi*180 << " " << f.r_s() << " "
	   << f.r_p() << " " << f.t_s() << " " << f.t_p() << " "
	   << f.R_s() << " " << f.R_p() << " " << f.T_s() << " "
	   << f.T_p()
	   << std::endl;
    }
  }
}

#endif
