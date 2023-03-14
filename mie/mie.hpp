#ifndef flick_mie
#define flick_mie

#include <complex>
#include <memory>
#include "../numeric/constants.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/function.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/legendre/legendre.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  class size_parameter {
    double sp_{1};
    const double pi_{constants::pi};
  public:
    double host_wavelength(double radius) {
      return 2 * pi_ * radius / sp_;
    }
    double radius(double host_wavelength) {
      return sp_ * host_wavelength / (2 * pi_);
    }
  };
  
  class basic_monodispeersed_mie {
  protected:
    std::complex<double> m_host_;
    std::complex<double> m_sphere_;
    double vacuum_wavelength_;
    const double pi_{constants::pi};
    double radius_{1e-6};
    stdvector angles_{0,pi_/2,pi_};
    const double host_wavelength_{m_host_.real()/vacuum_wavelength_};
    const std::complex<double> relative_refractive_index_{m_sphere_/m_host_};
    int precision_{2};
    double size_parameter() const {
      return 2 * pi_ * radius_ / host_wavelength_; 
    }
  public:
    basic_monodispeersed_mie(std::complex<double> m_host,
	      std::complex<double> m_sphere,
	      double vacuum_wavelength)
      : m_host_{m_host}, m_sphere_{m_sphere}, vacuum_wavelength_{vacuum_wavelength} {
    }
    virtual void radius(double r) = 0;
    virtual void angles(const stdvector& angles) = 0;
    virtual double extinction_cross_section() const = 0;
    virtual double absorption_cross_section() const = 0;
    virtual double asymmetry_factor() const = 0;
    virtual stdvector scattering_function(size_t row, size_t col) const = 0;
    void precision(size_t n) {
      precision_ = n;
    }
    int precision() const {
      return precision_;
    }
  };
   
  class parameterized_monodisperesed_mie : public basic_monodispeersed_mie
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
      basic_monodispeersed_mie::basic_monodispeersed_mie(m_host,
							 m_sphere,
							 vacuum_wavelength) {
      g0_ = read<pl_function>("mie/g_parameterized.txt");
    }
    void radius(double r) {
      radius_ = r;
      update_efficiency();
    }
    void angles(const stdvector& angles) {
      angles_ = angles;
    }
    double extinction_cross_section() const {
      return 2 * geometrical_cross_section();
    }
    double absorption_cross_section() const {
      return absorption_efficiency_ * geometrical_cross_section();
    }
    stdvector scattering_function(size_t row, size_t col) const {
      if (row == 0 && col == 0) {
	stdvector p(angles_.size());
	double g = asymmetry_factor();
	for (size_t i=0; i<p.size(); ++i) {
	  p[i] = henyey_greenstein(g).phase_function(angles_[i]);
	}
	return p;
      }
      else
	return stdvector(angles_.size(),0);
    }    
  };
  
  class size_distribution {
  protected:
    const double pi_{constants::pi};
    double a_;
    double b_;
  public:
    size_distribution(double a, double b)
      : a_{a}, b_{b} {} 
    virtual double center() const = 0;
    virtual double width() const = 0;
    virtual double value(double x) const = 0;
    virtual double weighted_integral(double alpha) const = 0;
  };
  
  class log_normal_distribution : public size_distribution
  // https://en.wikipedia.org/wiki/Log-normal_distribution
  {
  public:
    log_normal_distribution(double mu, double sigma)
      : size_distribution(mu,sigma) {
    }
    double value(double x) const {
      return 1/(x*b_*sqrt(2*pi_))*exp(-pow(log(x)-a_,2)/(2*pow(b_,2))); 
    }
    double center() const {
      return exp(a_-pow(b_,2));
    }
    double width() const {
      return b_;
    }
    double weighted_integral(double alpha) const {
      double k = 1/(2*pow(b_,2));
      double e = pow(2*k*a_+alpha,2)/(4*k) - k*pow(a_,2);
      std::cout << "k " << k<<std::endl;
      std::cout << "e " << e<<std::endl;
      //return 1/(b_*sqrt(2*pi_))*sqrt(pi_/k)*exp(e);
      double arg = pow(a_/b_ + alpha*b_,2);
      std::cout << "arg " << arg<<std::endl;
      return exp(0.5*(arg - pow(a_/b_,2)));
    }

  };
  /*
  class geometric_distribution : public size_distribution {
  };

  class gamma_distribution : public size_distribution {
  };
  */
  class basic_quantity {
  protected:
    basic_monodispeersed_mie& bm_;
    const size_distribution& sd_;
    stdvector center_quantity_;
    double alpha_{0};
  public:
    basic_quantity(basic_monodispeersed_mie& bm,
		   const size_distribution& sd)
      : bm_{bm}, sd_{sd} {}
    virtual size_t size() const = 0;
    virtual stdvector value(double x) = 0;
    stdvector center_quantity() const {
      return center_quantity_;
    }
    double distribution_center() const {
      return sd_.center();
    }
    double width() const {
      return sd_.width();
    }
    const size_distribution& sd() const {
      return sd_;
    }
    double alpha() {
      return alpha_;
    }
  };
  
  struct extinction_quantity : public basic_quantity {
    extinction_quantity(basic_monodispeersed_mie& bm,
			     const size_distribution& sd)
      : basic_quantity(bm,sd) {
      double rc = sd_.center();
      bm_.radius(rc);
      alpha_ = 2;
      center_quantity_ = {bm_.extinction_cross_section()/pow(rc,alpha_)};
    }
    size_t size() const {
      return 1;
    }
    stdvector value(double x) {
      double r = exp(x);
      bm_.radius(r);
      double rc = sd_.center();
      double v = bm_.extinction_cross_section()*r*sd_.value(r)
      - center_quantity_[0]*pow(r,alpha_)*r*sd_.value(r);
      return stdvector{v};
    }
  };
 
  class polydispersed_mie {
    std::shared_ptr<basic_quantity> bq_;
    basic_monodispeersed_mie& bm_;
    const size_distribution& sd_;
     size_t n_q_points_{0};
    double relative_rms_error(const stdvector& v1, const stdvector& v2) {
      return rms(v1/(v2+std::numeric_limits<double>::epsilon()));
    }
    stdvector integral(double x1, double x2,
		       const stdvector& compare_a) {
      std::vector<size_t> n_points{2,4,8,16,32,64,100};
      double error = std::numeric_limits<double>::max();
      double error_goal = pow(10,-bm_.precision());
      stdvector a(bq_->size(),0);
      stdvector previous_a(bq_->size(),0);
      size_t n = 0;
      while (error > error_goal && n < n_points.size()) {
	a = gl_integral_vector(*bq_,n_points[n]).value(x1,x2);
	n_q_points_ += n_points[n];
	if (n > 0) {
	  error = relative_rms_error(a-previous_a, compare_a);
	
	std::cout << "  compare_a "<< compare_a[0] << std::endl;
	std::cout << "  a "<< a[0] << std::endl;
	std::cout << "  da "<<  a-previous_a << std::endl;
	std::cout << "  previous a "<< previous_a[0] << std::endl;
	std::cout << "  error "<< error << std::endl;
	std::cout << "  Using "<<n_points[n] << " quadrature points\n";
	}
	previous_a = a;
	n++;	
      }
      if (n==n_points.size()) {
	std::cout << "Warning: Mie sub integral still at "<< error*100
		  << " % error using 100 quadrature points";
      }
      return a;
    }  
    stdvector integral() {   
      double error_goal = std::pow(10,-bm_.precision());
      double x0 = log(bq_->distribution_center());
      std::cout << "x0: "<< x0<< std::endl;
      std::cout << "r0: "<< exp(x0)<< std::endl;
      stdvector previous_a = bq_->center_quantity()
	* bq_->sd().weighted_integral(bq_->alpha());
      stdvector a = previous_a;
      std::cout << "previous a: " << previous_a[0]<< std::endl;
      stdvector direction{-1, 1};
      std::cout << "total area "<<a[0] << std::endl;
      for (size_t i=0; i<2; ++i) { // Both sides of max
	double error = std::numeric_limits<double>::max();
	std::cout << "Direction !!!!!!!!!!!!!!!!!!!!!!!! " << direction[i] << std::endl;
	double x1 = x0;
	size_t n = 0;
	while (error > error_goal) {
	  std::cout << "error: "<<error <<", goal "
		    << error_goal << std::endl;
	  double x2 = x1 + bq_->width()*direction[i];	 
	  std::cout << "x1 and x2: "<<x1 <<" " << x2<< std::endl;
	  a += integral(x1, x2, previous_a)*direction[i];
	  std::cout << "total area "<<a[0] << std::endl;
	  std::cout << "a: " << a[0]<< std::endl;
	  error = relative_rms_error(a-previous_a,a);
	  previous_a = a;
	  x1 = x2;
	  n++;
	}
      }
      std::cout << "Used a total of " << n_q_points_
      << " quadrature points"<<std::endl;
      return a;
    }
  public:
    polydispersed_mie(basic_monodispeersed_mie& bm,
		      const size_distribution& sd)
      : bm_{bm}, sd_{sd} {}
    double extinction_cross_section() {
      bq_ = std::make_shared<extinction_quantity>(bm_,sd_);
      return integral()[0];
    }
  };

  // class parameterized_polydisperesed_mie : public basic_polydisperesed_mie {
  //public:
  //   basic_basic_polydisperesed_mie()
  // };
 
  //};
  /*
  //template<class Size_distribution, class Monodispersed_mie>
  class polydispersed_mie {
    size_t precision_{2};
    const size_distribution& sd_;
    basic_monodispeersed_mie& mm_;
    //bool has_changed() {
    //}
    class basic_quantity {
    protected:
      stdvector q_;
      stdvector q_center_;
      basic_monodispeersed_mie& bm_;
      double x_center_;
    public:
      basic_quantity(basic_monodispeersed_mie& bm, size_t n, double x_center)
	: bm_{bm}, q_(n), q_center_(n), x_center_{x_center} {
	//q_center
      }
      stdvector& center_value() {
	return q_center_;
      }
    };
    class extinction : public basic_quantity {
      //basic_quantity(basic_monodispeersed_mie& bm, size_t n, double x_center) {
      using basic_quantity::basic_quantity;
      //q_center_[0] = mm_.extinction_coefficient(x_center_);
      //}
    };

    double approximate_integrangd(double x) {
      return 0;//center_value()[0];// * sd_.value(x);
    }
    double rest_integrand(double x) {
      return 0;
    }
    double integral(double x1, double x2) {
      return 0;
    }
  public:
    polydispersed_mie(const size_distribution& sd,
		      basic_monodispeersed_mie& mm)
      : sd_{sd}, mm_{mm} {}
    void set_precision(size_t n) {
      precision_ = n;
    }
    double extinction_cross_section() {
      return 0; //integral(extinction(mm_))[0];
    //  
    //  double a = 0;
    //  double r_1 = sd_.mean();
      //double r_2 = r_1 + sd_.width()*direction;
      */
      /*
      while (has_changed()) {
	weighted_iop_function f(mm_,sd_,"extinction_cross_section");
	da = integral(f,r1,r2);
	r_1 = r_2;
	r_2 = next_radius(r1, da);
	}
*/
  //  }
  //};
 
  
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

void expansion_coefficients(double m, double x1, int n_max, stdvector& a, stdvector& b) {
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
