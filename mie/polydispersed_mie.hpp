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
      return exp(a_+3*pow(b_,2)); // median size volume distribution
    }
    double width() const {
      return b_;
    }
    double weighted_integral(double alpha) const {
      return exp(a_*alpha+0.5*pow(alpha*b_,2));
    }
    double average_area() {
      return constants::pi*weighted_integral(2);
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
      return center_quantity_.size();
    }
    stdvector transformed_value(const stdvector& quantity) {
      double r = bm_.radius();
      return quantity * r * sd_.value(r);
      //return (quantity - center_quantity_ * pow(r,alpha_)) * r * sd_.value(r);
    }
    stdvector transformed_center(const stdvector& quantity) {
      return 0*quantity;
      //return quantity/pow(bm_.radius(),alpha_);
    }
  };
 
  struct absorption_quantity : public basic_quantity {
    absorption_quantity(basic_monodispersed_mie& bm,
			     const size_distribution& sd)
      : basic_quantity(bm,sd) {
      alpha_ = 3;
      bm_.radius(sd_.center());
      center_quantity_ = transformed_center({bm_.absorption_cross_section()});
    }  
    stdvector value(double x) {
      bm_.radius(exp(x));
      return transformed_value({bm_.absorption_cross_section()});
    }
  };

  struct scattering_quantity : public basic_quantity {
    scattering_quantity(basic_monodispersed_mie& bm,
			const size_distribution& sd)
      : basic_quantity(bm,sd) {
      alpha_ = 2;
      bm_.radius(sd_.center());
      center_quantity_ = transformed_center({bm_.scattering_cross_section()});
    }  
    stdvector value(double x) {
      bm_.radius(exp(x));
      return transformed_value({bm_.scattering_cross_section()});
    }
  };

  struct smatrix_quantity : public basic_quantity {
    size_t row_;
    size_t col_;
    smatrix_quantity(basic_monodispersed_mie& bm,
		     const size_distribution& sd,
		     size_t row, size_t col)
      : basic_quantity(bm,sd), row_{row}, col_{col} {
      alpha_ = 2;
      bm_.radius(sd_.center());
      center_quantity_ = transformed_center({bm_.scattering_matrix_element(row_,col_)});
    }  
    stdvector value(double x) {
      bm_.radius(exp(x));
      return transformed_value(bm_.scattering_matrix_element(row_,col_));
    }
  };
 
  template<class Monodispersed_mie, class Size_distribution>
  class polydispersed_mie {
    const double epsilon_ = std::numeric_limits<double>::epsilon()*10;
    std::shared_ptr<basic_quantity> bq_;
    Monodispersed_mie mm_;
    Size_distribution sd_;
    int precision_{2};
    
    size_t n_quadrature_points_;
    pl_function xy_points_;
    stdvector x_;
    stdvector y_;

    double error_goal() {
      return pow(10,-precision_+1);
    }
    double relative_area_change(const stdvector& da, const stdvector& a,
				double dx) {
      return vec::rms(da/dx*bq_->sd().width()/a);
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
	  error = relative_area_change(a-previous_a,compare_a+a,x2-x1);
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
      stdvector direction{-1, 1};
      xy_points_.clear();
      double max_step_width = 0.4;
      double step_width = max_step_width;
      for (size_t i=0; i<2; ++i) { // Both sides of max
	double error = std::numeric_limits<double>::max();
	double x1 = x0;
	while (error > error_goal()) {
	  double x2 = x1 + step_width * bq_->sd().width()*direction[i];	
	  stdvector da = integral(x1, x2, a)*direction[i];
	  total_q_points += n_quadrature_points_;
	  
	  if ( n_quadrature_points_ > 90) {
	    step_width /= 2;
	  } else {
	    a += da;

	    auto [x,y]=gl_integral_vector(*bq_, n_quadrature_points_).
	      xy_integration_points(x1,x2);	    
	    xy_points_ = concatenate(pl_function(x,y[0]),xy_points_);
	    
	    error = relative_area_change(da,a+da,x2-x1);
	    x1 = x2;
	  }	  
	  if (n_quadrature_points_ < 8) {
	    step_width *= 2;
	  }
	  if (step_width > max_step_width)
	    step_width = max_step_width;
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
    polydispersed_mie(Monodispersed_mie mm, Size_distribution sd)
      : mm_{mm}, sd_{sd} {
    }
    void precision(size_t n) {
      precision_ = n;
    }
    pl_function xy_points() {
      return xy_points_;
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
    double scattering_efficiency() {
      return scattering_cross_section()/sd_.average_area();
    }
    double absorption_efficiency() {
      return absorption_cross_section()/sd_.average_area();
    }

  };
}

#endif
