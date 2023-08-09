#ifndef flick_polydispersed_mie
#define flick_polydispersed_mie

#include "basic_monodispersed_mie.hpp"
#include "../numeric/legendre/legendre.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/physics_function.hpp"
#include "../environment/input_output.hpp"
#include <tuple>

namespace flick {
  class basic_quantity {
  protected:
    basic_monodispersed_mie& bm_;
    const size_distribution& sd_;
    stdvector center_quantity_;
    double alpha_{0};
    size_t size_{1};
    bool do_center_subtraction_ = true;
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
      if (do_center_subtraction_)
	return (quantity - center_quantity_ * pow(r,alpha_)) * r * sd_.value(r);
      return quantity * r * sd_.value(r);
    }
    stdvector transformed_center(const stdvector& quantity) {
      if (do_center_subtraction_)
	return quantity/pow(bm_.radius(),alpha_);
      return 0*quantity;
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
    const double epsilon_ = std::numeric_limits<double>::epsilon();
    std::shared_ptr<basic_quantity> bq_;
    Monodispersed_mie mm_;
    Size_distribution sd_;
    double accuracy_{0.05};
    stdvector integral_of_abs_integrand_;
    size_t n_quadrature_points_;
    size_t last_quadrature_element_ = 0;
    pl_function xy_points_;

    double relative_area_change(const stdvector& da, const stdvector& a,
				double dx) {
      double w = 5*bq_->sd().width();
      double change = vec::rms(sqrt(w/std::fabs(dx))*da/a);
      if (not isfinite(change))
      	return 0;           
      return change;
    }
    stdvector integral(double x1, double x2,
		       const stdvector& compare_a) {
      std::vector<size_t> n_points{1,2,4,8,16,32,64,128,256,512,1024,2048};
      double error = std::numeric_limits<double>::max();    
      stdvector a(bq_->size(),0);
      stdvector previous_a(bq_->size(),0);
      if (last_quadrature_element_ > 1)
	last_quadrature_element_--;
      size_t n = last_quadrature_element_;
      while (error > accuracy_ and n < n_points.size()) {
	size_t log2_n_points = log(n_points[n])/log(2);
	gl_integral_vector gl(bq_,log2_n_points);
	a = gl.value(x1,x2);
	n_quadrature_points_ = n_points[n];
	if (n > 0) {
	  error = relative_area_change(a-previous_a,compare_a+a,x2-x1);
	}
	previous_a = a;
	n++;
	integral_of_abs_integrand_ = gl.of_abs_integrand(x1,x2);
      }
      return a;
    }  
    stdvector integral() {
      double x0 = log(bq_->sd().center());
      stdvector a = bq_->center_quantity()
      	* bq_->sd().weighted_integral(bq_->alpha());
      stdvector direction{-1, 1};
      xy_points_.clear();
      double max_step_factor = 0.25;
      double step_factor = max_step_factor;
      for (size_t i=0; i<2; ++i) { // Both sides of max
	double error = std::numeric_limits<double>::max();
	double x1 = x0;
	while (error > accuracy_) {
	  double x2 = x1 + step_factor * bq_->sd().width()*direction[i];
	  stdvector da = integral(x1, x2, a)*direction[i];
	  if (n_quadrature_points_ >= 2048) {
	    step_factor /= 2;
	  } else {
	    a += da;
	    size_t log2_n_points = log(n_quadrature_points_)/log(2);
	    auto [x,y]=gl_integral_vector(bq_, log2_n_points).
	      xy_integration_points(x1,x2);
	    xy_points_ = concatenate(pl_function(x,y[0]),xy_points_);
	    double dx = x2-x1;
	    stdvector da_conservative = integral_of_abs_integrand_;
	    error = relative_area_change(da_conservative,a,dx);
	    x1 = x2;
	  }
	  if (n_quadrature_points_ < 8) {
	    step_factor *= 2;
	  }
	  if (step_factor > max_step_factor)
	    step_factor = max_step_factor;
	}
      }
      return a;
    }
  public:
    polydispersed_mie(Monodispersed_mie mm, Size_distribution sd)
      : mm_{mm}, sd_{sd} {
    }
    void percentage_accuracy(double p) {
      accuracy_ = p/100;
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
    stdvector normalized_scattering_matrix_element(size_t row, size_t col)
    /* Note that the normalization is such that the upper left corner
       element is normalized to 4*pi, not one, which is normally the
       flick convention for the phase function */
    {
      double k = 4*constants::pi/scattering_cross_section();
      return k*scattering_matrix_element(row,col);
    }
    std::tuple<std::vector<stdvector>, std::vector<stdvector>, stdvector> ab_functions() {
      std::vector<stdvector> a(4);
      a[0] = normalized_scattering_matrix_element(0,0);
      a[1] = normalized_scattering_matrix_element(1,1);
      a[2] = normalized_scattering_matrix_element(2,2);
      a[3] = normalized_scattering_matrix_element(3,3);
      std::vector<stdvector> b(2);
      b[0] = normalized_scattering_matrix_element(0,1);
      b[1] = normalized_scattering_matrix_element(2,3);
      stdvector x = vec::cos(mm_.angles());
      std::reverse(x.begin(),x.end());
      for (size_t i=0; i<a.size(); ++i)
	std::reverse(a[i].begin(),a[i].end());
      for (size_t i=0; i<b.size(); ++i)
	std::reverse(b[i].begin(),b[i].end());
      return {a,b,x};
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
