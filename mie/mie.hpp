#ifndef flick_mie
#define flick_mie

#include "../numeric/constants.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/function.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/legendre/legendre.hpp"
#include "../numeric/named_bounded_types.hpp"
#include "../environment/input_output.hpp"
#include <complex>

namespace flick {
  using stdcomplex = std::complex<double>;

  class refractive_index {
    stdcomplex value_{1,0};
  public:
    refractive_index() = default;
    refractive_index(const stdcomplex& m) : value_{m} {
      bounded_type<double, zero, std::exa>(m.real());
      bounded_type<double, zero, std::exa>(m.imag());
    }
    refractive_index(double real, double imag)
      : refractive_index{stdcomplex{real,imag}} {
    }
    const stdcomplex& value() const {
      return value_;
    }
    stdcomplex wavenumber(double vacuum_wavelength) const {
      return 2 * constants::pi * value_ / vacuum_wavelength;
    }
    stdcomplex size_parameter(double vacuum_wavelength, double radius) const {
      return wavenumber(vacuum_wavelength) * radius;
    }
    stdcomplex relative_to(const refractive_index& m) const {
      return value_ / m.value();
    }
  };
    
  class basic_monodispeersed_mie {
  protected:
    const double pi_{constants::pi};
    
    refractive_index m_host_;
    refractive_index m_sphere_;
    wavelength vacuum_wl_;
    double radius_{1e-6};
    stdvector angles_{0,pi_/2,pi_};

    int precision_{3};
  public:
    basic_monodispeersed_mie(const refractive_index& m_host,
			     const refractive_index& m_sphere,
			     const wavelength& vacuum_wl)
      : m_host_{m_host}, m_sphere_{m_sphere}, vacuum_wl_{vacuum_wl} {
    }
    virtual void radius(double r) = 0;
    virtual void angles(const stdvector& angles) = 0;
    virtual double absorption_cross_section() const = 0;
    virtual double scattering_cross_section() const = 0;
    virtual stdvector scattering_matrix_element(size_t row, size_t col) const = 0;

    void precision(size_t n) {
      precision_ = n;
    }
    int precision() const {
      return precision_;
    }
    double radius() const {
      return radius_;
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
    double Qa_{0};
    pl_function g0_;
    double n_ = m_sphere_.relative_to(m_host_).real();

    void update_efficiency() {    
      stdcomplex arg = 1./n_ * (pow(n_,3) - pow(pow(n_,2)-1., 3./2));
      double Qa0 = 8./3 * m_sphere_.value().imag()
	* m_host_.size_parameter(vacuum_wl_.value(),radius_).real()
	* std::abs(arg);
      Qa_  = 0.94 * (1 - exp(-Qa0 / 0.94));
    }
    double geometrical_cross_section() const {
      return pi_ * pow(radius_,2);
    }
    double asymmetry_factor() const {
      return pow(g0_.value(n_), pow(1-Qa_, 0.6));
    }
    double extinction_cross_section() const {
      return 2 * geometrical_cross_section();
    }
  public:
    parameterized_monodisperesed_mie(const refractive_index& m_host,
				     const refractive_index& m_sphere,
				     const wavelength& vacuum_wl) :
      basic_monodispeersed_mie::basic_monodispeersed_mie(m_host,
							 m_sphere,
							 vacuum_wl) {
      g0_ = read<pl_function>("mie/g_parameterized.txt");
    }
    void radius(double r) {
      radius_ = r;
      update_efficiency();
    }
    void angles(const stdvector& angles) {
      angles_ = angles;
    }
    double absorption_cross_section() const {
      return Qa_ * geometrical_cross_section();
    }
    double scattering_cross_section() const {
      return extinction_cross_section() - absorption_cross_section();
    }
    stdvector scattering_matrix_element(size_t row, size_t col) const
    // Note that integratinig element 0,0 over all 4*pi solid angles gives
    // the scattering cross section, where we count from 0 instead of one.
    {
      if (row == 0 && col == 0) {
	stdvector hg(angles_.size());
	double g = asymmetry_factor();
	for (size_t i=0; i<hg.size(); ++i) {
	  hg[i] = henyey_greenstein(g).phase_function(angles_[i]);
	}
	return hg * scattering_cross_section();
      }
      else
	return stdvector(angles_.size(),0);
    }    
  };

  class monodisperesed_mie : public basic_monodispeersed_mie
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
      double arg = pow(a_/b_ + alpha*b_,2);
      return exp(0.5*(arg - pow(a_/b_,2)));
    }
  };
  
  class geometric_distribution : public size_distribution {
  };

  class gamma_distribution : public size_distribution {
  };
  
  class basic_quantity {
  protected:
    basic_monodispeersed_mie& bm_;
    const size_distribution& sd_;
    stdvector center_quantity_;
    double alpha_{0};
    size_t size_{1};
  public:
    basic_quantity(basic_monodispeersed_mie& bm,
		   const size_distribution& sd)
      : bm_{bm}, sd_{sd} {
    }
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
    double alpha() const {
      return alpha_;
    }
    size_t size() const {
      return size_;
    }
  };
 
  struct absorption_quantity : public basic_quantity {
    absorption_quantity(basic_monodispeersed_mie& bm,
			     const size_distribution& sd)
      : basic_quantity(bm,sd) {
      double rc = sd_.center();
      bm_.radius(rc);
      alpha_ = 3;
      center_quantity_ = {bm_.absorption_cross_section()/pow(rc,alpha_)};
      size_ = 1;
    }  
    stdvector value(double x) {
      double r = exp(x);
      bm_.radius(r);
      double v = (bm_.absorption_cross_section() - center_quantity_[0]
		  * pow(r,alpha_)) * r * sd_.value(r);
      return stdvector{v};
    }
  };

  struct scattering_quantity : public basic_quantity {
    scattering_quantity(basic_monodispeersed_mie& bm,
			const size_distribution& sd)
      : basic_quantity(bm,sd) {
      double rc = sd_.center();
      alpha_ = 2;
      bm_.radius(rc);
      center_quantity_ = {bm_.scattering_cross_section()/pow(rc,alpha_)};
      size_ = 1;
    }  
    stdvector value(double x) {
      double r = exp(x);
      bm_.radius(r);
      double v = (bm_.scattering_cross_section() - center_quantity_[0]
	* pow(r,alpha_)) * r * sd_.value(r);
      return stdvector{v};
    }
  };

  struct smatrix_quantity : public basic_quantity {
    size_t row_;
    size_t col_;
    smatrix_quantity(basic_monodispeersed_mie& bm,
		     const size_distribution& sd,
		     size_t row, size_t col)
      : basic_quantity(bm,sd), row_{row}, col_{col} {
      double rc = sd_.center();
      alpha_ = 0;
      bm_.radius(rc);
      center_quantity_ = bm_.scattering_matrix_element(row_,col_)
	/ pow(rc,alpha_);
      size_ = center_quantity_.size();
    }  
    stdvector value(double x) {
      double r = exp(x);
      bm_.radius(r);
      return (bm_.scattering_matrix_element(row_,col_)
	      - center_quantity_*pow(r,alpha_)) * r * sd_.value(r);
    }
  };
  
  class polydispersed_mie {
    std::shared_ptr<basic_quantity> bq_;
    basic_monodispeersed_mie& bm_;
    const size_distribution& sd_;

    const double error_goal_ = pow(10,-bm_.precision());

    size_t n_quadrature_points_;

    double relative_rms_error(const stdvector& v1, const stdvector& v2) {
      return rms(v1/(v2+std::numeric_limits<double>::epsilon()));
    }
    stdvector integral(double x1, double x2,
		       const stdvector& compare_a) {
      std::vector<size_t> n_points{2,4,8,16,32,64,100};
      double error = std::numeric_limits<double>::max();
      
      stdvector a(bq_->size(),0);
      stdvector previous_a(bq_->size(),0);
      size_t n = 0;
      while (error > error_goal_ && n < n_points.size()) {
	a = gl_integral_vector(*bq_,n_points[n]).value(x1,x2);	
	n_quadrature_points_ = n_points[n];
	if (n > 0) {
	  error = relative_rms_error(a-previous_a, compare_a);

	  /*
	  std::cout << "  compare area "<< compare_a << std::endl;
	  std::cout << "  rest area "<< a << std::endl;
	  std::cout << "  delta area "<<  a-previous_a << std::endl;
	  std::cout << "  current error "<< error << std::endl;
	  std::cout << "  using "<<n_points[n] << " quadrature points\n";
	  */
	  
	}
	previous_a = a;
	n++;	
      }
      return a;
    }  
    stdvector integral() {
      size_t total_q_points = 0;
      size_t n_intervals = 0;
      double x0 = log(bq_->distribution_center());
      stdvector previous_a = bq_->center_quantity()
	* bq_->sd().weighted_integral(bq_->alpha());
      stdvector a = previous_a;
      stdvector direction{-1, 1};

      //std::cout << "\nEstimated total area before loop "<<a << std::endl;

      for (size_t i=0; i<2; ++i) { // Both sides of max
	double width_factor = 3;
	double error = std::numeric_limits<double>::max();

	//std::cout << "Direction " << direction[i] << std::endl;

	double x1 = x0;
	while (error > error_goal_) {
	  double x2 = x1 + width_factor * bq_->width()*direction[i];	 

	  //std::cout << "x1 and x2: "<<x1 <<" " << x2<< std::endl;

	  stdvector da = integral(x1, x2, previous_a)*direction[i];
	  total_q_points += n_quadrature_points_;	  
	  if (n_quadrature_points_ > 90) {
	    width_factor /= 2;
	  } else {
	    a += da;

	    //std::cout << "total area "<< a << std::endl;

	    error = relative_rms_error(a-previous_a,a);
	    previous_a = a;
	    x1 = x2;
	  }
	  /*
	  std::cout << "Width factor: "<< width_factor << std::endl;
	  std::cout << "error: "<<error <<", goal "
	  	    << error_goal_ << std::endl;
	  */
	  if (n_quadrature_points_ < 8) {
	    width_factor *= 2;
	  }
	  n_intervals++;
	}
      }
      /*
      std::cout << "Used a total of " << total_q_points
      		<< " quadrature points and "<<n_intervals
      		<< " sub intervals." <<std::endl;
      */
      return a;
    }
  public:
    polydispersed_mie(basic_monodispeersed_mie& bm,
		      const size_distribution& sd)
      : bm_{bm}, sd_{sd} {}
    double absorption_cross_section() {
      bq_ = std::make_shared<absorption_quantity>(bm_,sd_);
      return integral()[0];
    }
    double scattering_cross_section() {
      bq_ = std::make_shared<scattering_quantity>(bm_,sd_);
      return integral()[0];
    }
    stdvector scattering_matrix_element(size_t row, size_t col) {
      bq_ = std::make_shared<smatrix_quantity>(bm_,sd_,row,col);
      return integral(); 
    }
  };
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


  Move to material to reduce deps.
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

