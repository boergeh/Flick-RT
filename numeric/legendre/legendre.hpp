#ifndef flick_legendre
#define flick_legendre

#include <vector>
#include "../../environment/input_output.hpp"
#include "../function.hpp"
#include "../std_operators.hpp"
#include "../wigner/wigner_d.hpp"

namespace flick {
  two_columns read_quadrature(size_t log2_n_points) {
    std::string fname ="numeric/legendre/gl_quadrature/q_"+
      std::to_string((int)pow(2,log2_n_points))+".txt";
    return read<two_columns>(fname);     
  }
      
  template<class Function>
  class basic_gl_integral {
  protected:
    two_columns quadrature_;
    std::shared_ptr<Function> f_;
  public:
    basic_gl_integral(const std::shared_ptr<Function>& f, size_t log2_n_points)
      : f_{f} {
      quadrature_ = read_quadrature(log2_n_points);
    }
  };
    
  template<class Function>
  class gl_integral : public basic_gl_integral<Function> {
  public:
    gl_integral(const std::shared_ptr<Function>& f, size_t log2_n_points)
     : basic_gl_integral<Function>(f,log2_n_points) {
    }
    double value(double from, double to) const {
      double range = to - from;
      double v = 0;
      for (size_t n = 0; n < this->quadrature_.size(); ++n) {
	double x = this->quadrature_.column(0)[n];
	double w = this->quadrature_.column(1)[n];
	double local_x = from + (x+1)/2*range;
	v += w * this->f_->value(local_x);
      }
      return v/2*range;
    }
  };

  template<class Function>
  class gl_integral_vector : public basic_gl_integral<Function> {
  public:
     gl_integral_vector(const std::shared_ptr<Function>& f, size_t log2_n_points)
     : basic_gl_integral<Function>(f,log2_n_points) {
    }
    std::tuple<stdvector,std::vector<stdvector>>
    xy_integration_points(double from, double to) {
      stdvector x = from + (this->quadrature_.column(0)+1)/2*(to-from);
      std::vector<stdvector> y(this->f_->size(),stdvector(x.size()));
      for (size_t i = 0; i < x.size(); ++i) {
	stdvector func = this->f_->value(x[i]);
	for (size_t j = 0; j < func.size(); ++j) {
	   y[j][i] = func[j];
	}
      }
      return {x,y};
    }
    stdvector value(double from, double to) {
      auto [x,y] = xy_integration_points(from,to);
      stdvector v(y.size());
      for (size_t i = 0; i < v.size(); ++i) {
	v[i] = vec::sum(this->quadrature_.column(1)*y[i])/2*(to-from);
      }
      return v;   
    }
    stdvector of_abs_integrand(double from, double to) {
      auto [x,y] = xy_integration_points(from,to);
      stdvector v(y.size());
      for (size_t i = 0; i < v.size(); ++i) {
	v[i] = vec::sum(this->quadrature_.column(1)*vec::abs(y[i]))/2*(to-from);
      }
      return v;   
    }
  };

  template<class Function>
  class accumulated_integral_vector {
    std::shared_ptr<Function> f_;
    double percent_accuracy_;
    size_t log2_n_points_{0};
    stdvector total_;
    stdvector previous_total_;
    const size_t log2_max_{11};
    bool has_converged_{false};
    bool has_converged_in_one_iteration_{false};
    bool keep_integration_points_{false};
    pl_function integration_points_;
  public:
    accumulated_integral_vector(const std::shared_ptr<Function>& f, double percent_accuracy)
      : f_{f}, percent_accuracy_{percent_accuracy}, total_(f->size(),0), previous_total_(f->size(),0) {
    }
    void set_total(const stdvector& value) {
      total_ = value;
      previous_total_ = value;
    }
    stdvector previous_total() const {
      return previous_total_;
    }
    bool has_converged() const {
      return has_converged_;
    }
    bool significant_added_value() const {
      return !has_converged_in_one_iteration_;
    }
    void keep_integration_points(bool b) {
      keep_integration_points_ = b;
    }
    pl_function integration_points() const {
      return integration_points_;
    }
    void add_value(double x1, double x2, double estimated_total_width=0) {
      double error = std::numeric_limits<double>::max();
      stdvector previous_a = stdvector(f_->size(),0);
      stdvector a = stdvector(f_->size(),0);
      size_t n = log2_n_points_;
      double f = 1;
      if (estimated_total_width > 0)
      	f = sqrt(estimated_total_width/fabs(x2-x1));
      while (error > percent_accuracy_ and n <= log2_max_) {
	a = gl_integral_vector(f_, n).value(x1, x2);
	stdvector abs_a = gl_integral_vector(f_, n).of_abs_integrand(x1, x2);
	error = 100 * f * vec::rms(2*(abs_a-previous_a)/(abs_a+previous_a+total_));
	if (not std::isfinite(error))
	  error = 0;
	previous_a = abs_a;
	n++; 
      }
      n--;
      update_integration_points(x1,x2,n);
      update_convergence(error, n);
      update_likely_needed_points(n);
      if (has_converged_) {
	previous_total_ = total_;
      }
      total_ += a;
    }
    stdvector total() const {
      return total_;
    }
    stdvector total(double x1, double x2) {
      add_value(x1,x2);
      return total();
    }
    void reset_convergence() {
      has_converged_ = false;
      has_converged_in_one_iteration_ = false;
    }
  private:
    void update_convergence(double error, size_t n) {	
      if (error < percent_accuracy_)
	has_converged_ = true;
      else
	has_converged_ = false;
      
      if (has_converged_ and n == log2_n_points_)
	has_converged_in_one_iteration_ = true;
      else
	has_converged_in_one_iteration_ = false;
    }
    void update_likely_needed_points(size_t n) {
      if (n > 2)
	log2_n_points_ = n-2;
    }
    void update_integration_points(double x1, double x2, size_t n) {
      if (keep_integration_points_ && has_converged_) {
	auto [x,y]=gl_integral_vector(f_, n).xy_integration_points(x1,x2);
	integration_points_ = concatenate(pl_function(x,y[0]),integration_points_);
      }
    }
  };

  template<class Function>
  class vector_function {
    std::shared_ptr<Function> f_;
  public:
    vector_function(const std::shared_ptr<Function>& f) : f_{f} {}
    size_t size() const {
      return 1;
    }
    stdvector value(double x) const {
      return stdvector{f_->value(x)};
    }
  };
  
  template<class Function>
  class accumulated_integral : public accumulated_integral_vector<vector_function<Function>> {
  public:
    accumulated_integral(const std::shared_ptr<Function>& f, double percent_accuracy) :
      accumulated_integral_vector<vector_function<Function>>(std::make_shared<vector_function<Function>>(f),percent_accuracy) {
    }
    void set_total(double value) {
      accumulated_integral_vector<vector_function<Function>>::set_total(stdvector{value});
    }
    double previous_total() const {
      return accumulated_integral_vector<vector_function<Function>>::previous_total()[0];
    }
    double total() const {
      return accumulated_integral_vector<vector_function<Function>>::total()[0];
    }
    double total(double x1, double x2) {
      return accumulated_integral_vector<vector_function<Function>>::total(x1,x2)[0];
    }
  };
  
  class legendre {
    std::vector<stdvector> p_;
  public:
    legendre(size_t n_terms, stdvector x) {
      p_.resize(x.size());
      for (size_t i=0; i < x.size(); ++i) {
	p_[i] = wigner_d(x[i],0,0,n_terms).terms();
      }
    }
    double value(size_t term_number, size_t x_number) {
      return p_[x_number][term_number];
    }
  };

  template<class Function>
  stdvector legendre_expansion(const std::shared_ptr<Function>& f,
			       size_t n_terms,
			       size_t log2_n_points = 4) {
    stdvector terms(n_terms);
    stdvector x = read_quadrature(log2_n_points).column(0);
    legendre legendre(n_terms, x);
    for (size_t i=0; i < terms.size(); ++i) {
      auto plf = std::make_shared<pl_function>();
      for (size_t j=0; j < x.size(); ++j) {
	double y = f->value(x[j]) * legendre.value(i,j);
	plf->append(point{x[j],y});
      }
      terms[i] = (2.*i+1)/2 * gl_integral(plf,log2_n_points).value(-1,1);
    }
    return terms;
  }

  class legendre_evaluation {
    const stdvector& coefficients_;
  public:
    legendre_evaluation(const stdvector& coefficients)
      : coefficients_{coefficients} {}
    stdvector values(const stdvector& x) {
      stdvector values(x.size());
      legendre legendre(coefficients_.size(), x);
      for (size_t i=0; i < x.size(); ++i) {
	double v = 0;
	for (size_t l=0; l < coefficients_.size(); ++l) {
	  v += coefficients_[l]*legendre.value(l,i);
	}
	values[i] = v;
      }
      return values;
    }
  };
}

#endif
