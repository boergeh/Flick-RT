#ifndef flick_polydispersed_mie
#define flick_polydispersed_mie

#include "basic_monodispersed_mie.hpp"
#include "../numeric/legendre/legendre.hpp"
#include "../numeric/std_operators.hpp"
#include "../environment/input_output.hpp"
#include <tuple>

namespace flick {
  class size_distribution
  // Size number distribution. Integral over all sizes equals one.
  {
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
    double particles_per_volume(double volume_fraction) const {
      return volume_fraction * 1/(4./3*pi_*weighted_integral(3));
    }
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
      return exp(a_*alpha+0.5*pow(alpha*b_,2));
    }
    static std::tuple<double, double>
    from_volume_distribution(double mu, double sigma) {
      return {mu-3*pow(sigma,2), sigma};
    }
  };
  
  //class geometric_distribution : public size_distribution {
  //};

  //class gamma_distribution : public size_distribution {
  //};
  
  class basic_quantity {
  protected:
    basic_monodispersed_mie& bm_;
    const size_distribution& sd_;
    stdvector center_quantity_;
    double alpha_{0};
    size_t size_{1};
  public:
    basic_quantity(basic_monodispersed_mie& bm,
		   const size_distribution& sd)
      : bm_{bm}, sd_{sd} {
    }
    virtual stdvector value(double x) = 0;

    stdvector center_quantity() const {
      return center_quantity_;
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
    absorption_quantity(basic_monodispersed_mie& bm,
			     const size_distribution& sd)
      : basic_quantity(bm,sd) {
      double rc = sd_.center();
      bm_.radius(rc);
      alpha_ = 2;
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
    scattering_quantity(basic_monodispersed_mie& bm,
			const size_distribution& sd)
      : basic_quantity(bm,sd) {
      double rc = sd_.center();
      alpha_ = 3;
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
    smatrix_quantity(basic_monodispersed_mie& bm,
		     const size_distribution& sd,
		     size_t row, size_t col)
      : basic_quantity(bm,sd), row_{row}, col_{col} {
      double rc = sd_.center();
      alpha_ = 1;
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
 
  template<class Monodispersed_mie, class Size_distribution>
  class polydispersed_mie {
    const double epsilon_ = std::numeric_limits<double>::epsilon()*10;
    std::shared_ptr<basic_quantity> bq_;
    Monodispersed_mie mm_;
    Size_distribution sd_;
    int precision_{3};
    size_t n_quadrature_points_;

    double error_goal() {
      return pow(10,-precision_+1);
    }
    double relative_rms_error(const stdvector& v1, const stdvector& v2) {
      return vec::rms(v1/(v2+std::numeric_limits<double>::epsilon()));
    }
    stdvector integral(double x1, double x2,
		       const stdvector& compare_a) {
      std::vector<size_t> n_points{2,4,8,16,32,64,100};
      double error = std::numeric_limits<double>::max();
      
      stdvector a(bq_->size(),0);
      stdvector previous_a(bq_->size(),0);
      size_t n = 0;
      while (error > error_goal() && n < n_points.size()) {
	a = gl_integral_vector(*bq_,n_points[n]).value(x1,x2);	
	n_quadrature_points_ = n_points[n];
	if (n > 0) {
	  error = relative_rms_error(a-previous_a, compare_a);
	  
	  //std::cout << "  compare area "<< compare_a << std::endl;
	  //std::cout << "  rest area "<< a << std::endl;
	  //std::cout << "  delta area "<<  a-previous_a << std::endl;
	  //std::cout << "  current error "<< error << std::endl;
	  //std::cout << "  using "<<n_points[n] << " quadrature points\n";	  
	}
	previous_a = a;
	n++;	
      }
      return a;
    }  
    stdvector integral() {
      size_t total_q_points = 0;
      size_t n_intervals = 0;
      double x0 = log(bq_->sd().center());
      stdvector a = bq_->center_quantity()
	* bq_->sd().weighted_integral(bq_->alpha());
      stdvector previous_a = a;
      stdvector direction{-1, 1};

      //std::cout << "\nEstimated total area before loop "<<a << std::endl;

      for (size_t i=0; i<2; ++i) { // Both sides of max
	double width_factor = 0.5;
	double error = std::numeric_limits<double>::max();

	//std::cout << "Direction " << direction[i] << std::endl;

	double x1 = x0;
	while (error > error_goal()) {
	  double x2 = x1 + width_factor * bq_->sd().width()*direction[i];	

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
	    //std::cout << "Width factor: "<< width_factor << std::endl;

	  }
	  
	  //std::cout << "error: "<<error <<", goal "
	  //	    << error_goal_ << std::endl;
	  
	  if (n_quadrature_points_ < 8) {
	    width_factor *= 2;
	  }
	  n_intervals++;
	}
      }
      
      //std::cout << "Used a total of " << total_q_points
      //		<< " quadrature points and "<<n_intervals
      //		<< " sub intervals." <<std::endl;
      
      return a;
    }
  public:
    polydispersed_mie(Monodispersed_mie mm, Size_distribution sd)
      : mm_{mm}, sd_{sd} {
    }
    void precision(size_t n) {
      precision_ = n;
    }
    double absorption_cross_section() {
      if (sd_.width() < epsilon_) {
	mm_.radius(sd_.center());
	return mm_.absorption_cross_section();
      } else {
	bq_ = std::make_shared<absorption_quantity>(mm_,sd_);
	return integral()[0];
      }
    }
    double scattering_cross_section() {
      if (sd_.width() < epsilon_) {
	mm_.radius(sd_.center());
	return mm_.scattering_cross_section();
      } else {
	bq_ = std::make_shared<scattering_quantity>(mm_,sd_);
	return integral()[0];
      }
    }
    stdvector scattering_matrix_element(size_t row, size_t col) {
      if (sd_.width() < epsilon_) {
	mm_.radius(sd_.center());
	return mm_.scattering_matrix_element(row,col);
      } else {
	bq_ = std::make_shared<smatrix_quantity>(mm_,sd_,row,col);
	return integral(); 
      }
    }
  };
}

#endif
