#ifndef flick_mie
#define flick_mie

#include <complex>
#include "../numeric/constants.hpp"
#include "../numeric/function.hpp"
#include "../numeric/physics_function.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  class basic_mie {
  protected:
    std::complex<double> m_host_;
    std::complex<double> m_sphere_;
    double vacuum_wavelength_;
    const double pi_{constants::pi};
    double radius_{1e-6};
    std::vector<double> angles_{0,pi_/2,pi_};
    const double host_wavelength_{m_host_.real()/vacuum_wavelength_};
    const std::complex<double> relative_refractive_index_{m_sphere_/m_host_};
    double size_parameter() const {
      return 2 * pi_ * radius_ / host_wavelength_; 
    }    
  public:
    basic_mie(std::complex<double> m_host,
	      std::complex<double> m_sphere,
	      double vacuum_wavelength)
      : m_host_{m_host}, m_sphere_{m_sphere}, vacuum_wavelength_{vacuum_wavelength} {
    }
    virtual void set_radius(double r) = 0;
    virtual void set_angles(const std::vector<double>& angles) = 0;
    virtual double extinction_cross_section() const = 0;
    virtual double absorption_cross_section() const = 0;
    virtual double asymmetry_factor() const = 0;
    virtual std::vector<double> scattering_function(size_t row, size_t col) const = 0;
  };
   
  class parameterized_monodisperesed_mie : public basic_mie
  // Approximate Mie-code solutions for large spheres, see:
  // Stamnes, K., Hamre, B., Stamnes, J.J., Ryzhikov, G., Biryulina,
  // M., Mahoney, R., Hauss, B. and Sei, A., 2011. Modeling of
  // radiation transport in coupled atmosphere-snow-ice-ocean
  // systems. Journal of Quantitative Spectroscopy and Radiative
  // Transfer, 112(4), pp.714-726.
  {
    double absorption_efficiency_{0};
    pl_function g0_;
    void update_efficiency() {
      double n = relative_refractive_index_.real();
      std::complex<double> arg = 1./n * (pow(n,3) - pow(pow(n,2)-1., 3./2));
      double Qa0=8./3*m_sphere_.imag()*size_parameter()*std::abs(arg);
      absorption_efficiency_  = 0.94 * (1 - exp(-Qa0 / 0.94));
    }
    double geometrical_cross_section() const {
      return pi_ * pow(radius_,2);
    }
    double asymmetry_factor() const {
      double n = relative_refractive_index_.real();
      double Qa = absorption_efficiency_;
      return pow(g0_.value(n), pow(1-Qa, 0.6));
    }
  public:
    parameterized_monodisperesed_mie(std::complex<double> m_host,
		      std::complex<double> m_sphere,
		      double vacuum_wavelength) :
      basic_mie::basic_mie(m_host,m_sphere,vacuum_wavelength) {
      g0_ = read<pl_function>("mie/g_parameterized.txt");
    }
    void set_radius(double r) {
      radius_ = r;
      update_efficiency();
    }
    void set_angles(const std::vector<double>& angles) {
      angles_ = angles;
    }
    double extinction_cross_section() const {
      double Qa = absorption_efficiency_;
      double Qs = 2 - Qa;
      return (Qa + Qs) * geometrical_cross_section();
    }
    double absorption_cross_section() const {
      return absorption_efficiency_ * geometrical_cross_section();
    }
    std::vector<double> scattering_function(size_t row, size_t col) const {
      if (row == 0 && col == 0) {
	std::vector<double> p(angles_.size());
	double g = asymmetry_factor();
	for (size_t i=0; i<p.size(); ++i) {
	  p[i] = henyey_greenstein(g).phase_function(angles_[i]);
	}
	return p;
      }
      else
	return std::vector<double>(angles_.size(),0);
    }
    
  };
  
  class size_distribution {
  public:
    virtual double mean() {
      return 0;
    }
    virtual double width() {
      return 0;
    }
  };
  class log_normal_distribution : public size_distribution {
    
  };

  class geometric_distribution : public size_distribution {
  };

  class gamma_distribution : public size_distribution {
  };

  template<class Size_distribution, class Monodisperesed_mie>
  class mie {
    size_t accuracy_;
    bool has_changed() {
    }
  public:
    mie(size_t accuracy) : accuracy_{accuracy} {}
    double extinction_cross_section() {
      //integrate(md,"")
      double a = 0;
      double r_1 = sd_.mean();
      double r_2 = r_1 + sd_.width()*direction;
      while (has_changed()) {
	weighted_iop_function f(mm_,sd_,"extinction_cross_section");
	da = integral(f,r1,r2);
	r_1 = r_2;
	r_2 = next_radius(r1, da);
      }
    }
  };
 
  
  /* Move to material to reduce deps.
  template<class Material>
  class material_refracitve_index : public refractive_index {
    const Material& m_;
  public:
    material_refractive_index(const Material& m) : m_{m} {}
    std::complex<double> value() {
      return std::complex<double>(m_.real_refractvie_index(),m_.imag...);
    } 
  };
  */

  /*
  class extinction_cross_section {
  public:
    value()
  }
  namespace mie {
    template<class Size_distribution, class Refractive_index, class Solver>
    class base {
    public:
      double extinction_cross_section() {
	return 0;
      }
      double scattering_cross_section() {
	return 0;
      }
      double scattering_function(size_t row, size_t col, double angle) {
	return 0;
      }
      void set_precision(size_t significant_figures) {
      }
    private:
      
      
    };
  }
  */
}
#endif

/*
#include <cmath>
#include <vector>

void expansion_coefficients(double m, double x1, int n_max, std::vector<double>& a, std::vector<double>& b) {
    a.resize(n_max);
    b.resize(n_max);

    for (int n = 1; n <= n_max; ++n) {
        double jmx = std::cyl_bessel_j(n+0.5, m*x1);
        double jx = std::cyl_bessel_j(n+0.5, x1);
        double hx = std::cyl_neumann(n+0.5, x1);
        double dmxjmx = m*x1*std::cyl_bessel_j(n-1+0.5, m*x1) - n*std::cyl_bessel_j(n+0.5, m*x1);
        double dxjx = x1*std::cyl_bessel_j(n-1+0.5, x1) - n*std::cyl_bessel_j(n+0.5, x1);
        double dxhx = x1*std::cyl_neumann(n-1+0.5, x1) - n*std::cyl_neumann(n+0.5, x1);  
        double A = jmx*dxjx;
        double B = jx*dmxjmx;
        double C = jmx*dxhx;
        double D = hx*dmxjmx;
        a[n-1] = (m*m*A-B)/(m*m*C-D);
        b[n-1] = (A-B)/(C-D);
    }
}

*/
