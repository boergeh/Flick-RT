#ifndef flick_function
#define flick_function
 
#include <optional>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <memory>
#include <stdexcept>
#include "sorted_vector.hpp"
#include "range.hpp"

namespace flick {
  using stdvec = std::vector<double>;
  
  std::string read_header(std::istream& is) {
    std::string h;
    std::string s;
    std::streampos start = is.tellg();
    is >> s;
    if (s.size() > 1 && s.substr(0,2)=="/*") {
      while (!(is.eof() || (s.size() > 1 && s.substr(s.size()-2)=="*/"))) {
	is >> s;
      }
      std::streamoff sz = is.tellg() - start; 
      is.seekg(start);
      h.resize(sz);
      is.read(&h[0], sz);
    } else
      is.seekg(start);
    return h+"\n\n";
  }

  struct point {
    double x_;
    double y_;
  public:
    point() = default;
    point(double x, double y) : x_(x), y_(y) {}
    double x() const {return x_;}
    double y() const {return y_;}
  };
 
  class basic_interpolation {
  protected:
    double x1, x2, y1, y2;
    const double epsilon = std::numeric_limits<double>::epsilon();
  public:
    basic_interpolation(const point& low, const point& high)
      : x1{low.x()}, x2{high.x()}, y1{low.y()}, y2{high.y()} {
    }
  };
  
  class piecewise_linear : public basic_interpolation {
    double a;
    double b;
  public:
    static step_type get_step_type() {
      return step_type::linear;
    }
    piecewise_linear(const point& low, const point& high)
      : basic_interpolation(low, high) {
      a = (y2-y1)/(x2-x1);
      b = y1-a*x1;
    }
    double y(double x) {
      return a * x + b;
    }
    double derivative(double x) {
      return a;
    }
    double integral(double limit_a, double limit_b) {
      return 0.5*a*(pow(limit_b,2)-pow(limit_a,2)) + b*(limit_b-limit_a); 
    }
    std::optional<double> integral_limit_b(double limit_a,
					   double integral_value) {
      if (fabs(a) < epsilon) {
	if (fabs(b) < epsilon)
	  return std::nullopt;
	return limit_a + integral_value/b;
      }
      double k = integral_value+0.5*a*pow(limit_a,2)+b*limit_a;
      double arg = pow(b,2)+2*a*k;
      if (arg < 0)
	return std::nullopt;
      if (integral_value > 0)
	return (-b+sqrt(arg))/a;
      return (-b-sqrt(arg))/a;
    }
  };

  class piecewise_exponential : public basic_interpolation {
    double a;
    double k;
  public:
    static step_type get_step_type() {
      return step_type::linear;
    }
    piecewise_exponential(const point& low, const point& high)
      : basic_interpolation(low, high) {
      k = log(y2/y1) / (x2-x1);
    }
    double y(double x) {
      if (y1 <= 0)
	return 0;
      return y1*pow(y2/y1,(x-x1)/(x2-x1));
    }
    double derivative(double x) {
      ensure(y1>0);
      return y(x)*k;
    }
    double integral(double limit_a, double limit_b) {
      if (y1 <= 0) 
	return 0;
      if (fabs(k) < epsilon)
	return y1*(limit_b-limit_a);
      return (y(limit_b)-y(limit_a))/k;
    }
    std::optional<double> integral_limit_b(double limit_a,
					   double integral_value) {
      if (y1 <= 0)
	return std::nullopt;
      if (fabs(k) < epsilon)
	return limit_a + integral_value / y1;
      double a = y1*exp(-k*x1);
      double arg = integral_value*k/a + exp(k*limit_a);
      if (arg <= 0)
	return std::nullopt;
      return log(arg)/k;
    }
  private:
    void ensure(bool b) {
      if (!b) {
	throw std::runtime_error("numeric piecewise exponential interpolation");
      }
    }    
  };

  class piecewise_power : public basic_interpolation {
    double k;
  public:
    static step_type get_step_type() {
      return step_type::exponential;
    }
    piecewise_power(const point& low, const point& high)
      : basic_interpolation(low, high) {
      k = log(y2/y1) / log(x2/x1);
    }
    double y(double x) {
      if (y1 <= 0)
	return 0;
      return y1*pow(x/x1,k);
    }
    double derivative(double x) {
      ensure(y1>0);
      if (fabs(k) < epsilon)
	return 0;
      return y1*k*pow(x/x1,k-1)/x1;
    }
    double integral(double limit_a, double limit_b) {
      if (y1 <=0 or limit_a <= 0 or limit_b <= 0)
	return 0;
      if (fabs(k+1) < epsilon)
	return y1*x1*log(limit_b/limit_a);
      return y1/(k+1)*x1*(pow(limit_b/x1,k+1)-pow(limit_a/x1,k+1));
    }
    std::optional<double> integral_limit_b(double limit_a,
					   double integral_value) {
      if (y1 <=0 or limit_a <= 0)
	return std::nullopt;
      if (fabs(k+1) < epsilon)
	return limit_a * exp(integral_value/(y1*x1));
      double arg = integral_value*(k+1)/(x1*y1) + pow(limit_a/x1,k+1);
      if (arg <= 0)
	return std::nullopt;
      return x1*pow(arg,1/(k+1));
    }
  private:
    void ensure(bool b) {
      if (!b) {
	throw std::invalid_argument("piecewise_power");
      }
    }    
  };
 
  template<class I>
  class function {
    std::string header_;
    sorted_vector xv_;
    stdvec yv_;
  public:
    sorted_vector xv2_;
    function() = default;
    function(double value) : xv_{{1}}, yv_{stdvec{value}} {}
    function(const stdvec& xv, const stdvec& yv)
    {
      if (yv.size()>1 && xv.front() > xv.back()) {
	stdvec a = xv;
	stdvec b = yv;
	std::reverse(a.begin(), a.end());
	std::reverse(b.begin(), b.end());
	xv_ = a;
	yv_ = b;
      } else {
	xv_ = xv;
	yv_ = yv;
      }	
      if (yv.size()>1) {
	ensure(xv.size()==yv.size());
	xv_.set_step_type(I::get_step_type());
      }
    }
    auto& clear() {
      xv_.clear();
      yv_.clear();
      return *this;
    }
    auto& add_extrapolation_points(double weight=1,double dx_scaling_factor=1) {
      ensure(xv_.size() > 1);
      if (I::get_step_type()!=step_type::linear)
	ensure(weight > 0);
      stdvec xv = xv_.all_values();
      double dx_front = (xv[1]-xv[0])*dx_scaling_factor;
      double dx_back = (xv[xv.size()-1]-xv[xv.size()-2])*dx_scaling_factor;
      if (I::get_step_type()!=step_type::linear) {
	dx_front = xv[0]*dx_front/xv[1];
	dx_back = xv.back()*dx_back/xv.back();
      }
      xv.insert(xv.begin(),xv[0]-dx_front/2);
      xv.insert(xv.begin(),xv[0]-dx_front);
      xv.emplace_back(xv.back()+dx_back/2);
      xv.emplace_back(xv.back()+dx_back);
      xv_ = sorted_vector{xv};
      yv_.insert(yv_.begin(),yv_[0]*weight);
      yv_.insert(yv_.begin(),yv_[1]*weight);
      yv_.emplace_back(yv_.back()*weight);
      yv_.emplace_back((yv_.end()[-2])*weight);
      return *this;
    }
    size_t size() const {
      return yv_.size();
    }
    auto& header(const std::string& h) {
      header_ = h;
      return *this;
    }
    std::string header() const {
      return header_;
    }
    auto& append(const point& p) {
      xv_.append(p.x());
      yv_.emplace_back(p.y());
      return *this;
    }
    auto& append(const stdvec& xv, const stdvec& yv) {
      ensure(xv.size()==yv.size());
      for (size_t i=0; i<xv.size(); ++i) {
	append(point{xv[i],yv[i]});
      }
      return *this;
    }
    auto& scale_x(double factor) {
      ensure(factor > 0);
      xv_.scale(factor);
      return *this;
    }   
    auto& scale_y(double factor) {
      ensure(factor > 0);
      for (size_t i=0; i<yv_.size(); ++i)
	yv_[i] *= factor;
      return *this;
    }
    auto& normalize() {
      scale_y(1/integral());
      return *this;
    }
    const stdvec& x() const {
      return xv_.all_values();
    }
    const stdvec& y() const {
      return yv_;
    }
    double value(double x=1) const {
      if (yv_.size()==1)
	return yv_[0];
      ensure(yv_.size() > 1);
      auto [p1, p2] = points_at(x);
      return I{p1,p2}.y(x);      
    }
    size_t low_index_near(double x) const {
      return xv_.find(x);
    }
    double derivative(double x) const {
      ensure(yv_.size() > 1);
      auto [p1, p2] = points_at(x);
      return I{p1,p2}.derivative(x);
    }
    std::optional<double> integral_limit_b(double limit_a, double integral_value) const {
      ensure(yv_.size() > 1);
      std::shared_ptr<sorted_vector::iterator> it;
      if (integral_value > 0)
	it = std::make_shared<sorted_vector::ascending_iterator>(xv_);
      else
	it = std::make_shared<sorted_vector::descending_iterator>(xv_);
      it->move_to_bin_at(limit_a);
      double area = 0;
      point p1;
      point p2;
      while(true) {
	p1 = previous_point(it.get());
	p2 = next_point(it.get());
	double next_area = I{p1,p2}.integral(limit_a, p2.x());
	if(fabs(area + next_area) > fabs(integral_value) || it->is_in_end_bin())
	  return I{p1,p2}.integral_limit_b(limit_a, integral_value-area);
	area += next_area;
	limit_a = p2.x();
	it->move_to_next_bin();
      }      
    }
    double integral(double limit_a, double limit_b) const {
      if (xv_.size()==0)
	return 0;
      if (xv_.size()==1)
	return yv_[0]*(limit_b-limit_a);
      std::shared_ptr<sorted_vector::iterator> it;
      if (limit_a < limit_b) {
	it = std::make_shared<sorted_vector::ascending_iterator>(xv_);
      } else {
	it = std::make_shared<sorted_vector::descending_iterator>(xv_);
      }
      it->move_to_bin_at(limit_a);
      double area = 0;
      point p1;
      point p2;
      while(true) {
	p1 = previous_point(it.get());
	p2 = next_point(it.get());
	if(p2.x() > limit_b || it->is_in_end_bin())
	  return area + I{p1,p2}.integral(limit_a, limit_b);
	area += I{p1,p2}.integral(limit_a, p2.x());
	limit_a = p2.x();
	it->move_to_next_bin();
      }      
    }
    double integral() const {
      return integral(xv_[0], xv_[xv_.size()-1]);
    }
    stdvec accumulation() const {
      ensure(xv_.size() > 1);
      stdvec a(xv_.size());
      a[0] = 0;
      double area = 0;
      for (size_t i=0; i < xv_.size()-1; ++i) {
	point p1 = {xv_[i], yv_[i]};
	point p2 = {xv_[i+1], yv_[i+1]};
	double da = I{p1,p2}.integral(p1.x(), p2.x());
	ensure(da >= 0);
	area += da;
	a[i+1] = area;
      }
      return a;
    }
  private:
    friend std::ostream& operator<<(std::ostream &os,
				    const function<I>& f) {
      os << f.header();
      for (size_t i = 0; i<f.yv_.size(); ++i)
	os << f.xv_[i] << " " << f.yv_[i] << '\n';
      return os;
    }
    friend std::istream& operator>>(std::istream &is,
				    function<I>& f) {
      f.header(read_header(is));
      double x, y;
      while(is >> x >> y) {
	f.append(point{x,y});
      }
      return is;
    }
    point next_point(sorted_vector::iterator *it) const {
      return point{xv_[it->next_index()],yv_[it->next_index()]};
    }
    point previous_point(sorted_vector::iterator *it) const {
      return point{xv_[it->previous_index()],yv_[it->previous_index()]};
    }
    std::tuple<point, point> points_at(double x) const {
      size_t n = xv_.find(x);
      point p1{xv_[n],yv_[n]};
      point p2{xv_[n+1],yv_[n+1]};
      return {p1,p2};
    }
    void ensure(bool b) const {
      if (!b)
	throw std::runtime_error("numeric function");
    }
  };

  template<class I>
  function<I> concatenate(function<I> fa, function<I> fb) {
    if (fa.x().empty())
      return fb;
    if (fb.x().empty())
      return fa;
    if (fa.x().back() < fb.x().front())      
      return fa.append(fb.x(),fb.y());
    else
      return fb.append(fa.x(),fa.y());
  }

  template<class I>
  function<I> add(const function<I>& fa, const function<I>& fb, const stdvec& xv) {
    stdvec yv(xv.size());
    for (size_t i=0; i < xv.size(); ++i)
      yv[i] = fa.value(xv[i])+fb.value(xv[i]);
    return function<I>{xv,yv};
  }

  template<class I>
  function<I> subtract(const function<I>& fa, const function<I>& fb, const stdvec& xv) {
    stdvec yv(xv.size());
    for (size_t i=0; i < xv.size(); ++i)
      yv[i] = fa.value(xv[i])-fb.value(xv[i]);
    return function<I>{xv,yv};
  }

  template<class I>
  function<I> scale_to_integral(function<I> f, double integral_value) {
    double fi = f.integral();
    if (fabs(fi-integral_value) > std::numeric_limits<double>::epsilon())
      return f.scale_y(integral_value/fi);
    return f;
  }

  template<class I>
  function<I> integral_conservative_add(const function<I>& fa, const function<I>& fb, const stdvec& xv) {
    function<I> f = add(fa,fb,xv);
    return scale_to_integral(f, fa.integral()+fb.integral());
  }

  template<class I>
  function<I> multiply(const function<I>& fa, const function<I>& fb, const stdvec& xv) {
    stdvec yv(xv.size());
    for (size_t i=0; i < xv.size(); ++i)
      yv[i] = fa.value(xv[i])*fb.value(xv[i]);
    return function<I>{xv,yv};
  }

  template<class I>
  function<I> divide(const function<I>& fa, const function<I>& fb, const stdvec& xv) {
    stdvec yv(xv.size());
    for (size_t i=0; i < xv.size(); ++i)
      yv[i] = fa.value(xv[i])/fb.value(xv[i]);
    return function<I>{xv,yv};
  }

  template<class I>
  function<piecewise_linear> inverted_cumulative_distribution(function<I> f) {
    f.normalize();
    return function<piecewise_linear>{f.accumulation(),f.x()};
  }
  
  template<class I>
  function<I> importance_sampled(const function<I>& f, size_t n_points) {
    function<piecewise_linear> inv_cum = inverted_cumulative_distribution(f);
    stdvec unit_interval = range(0,1,n_points).linspace();
    stdvec x(n_points);
    stdvec y(n_points);
    for (size_t i=0; i < n_points; ++i) {
      x[i] = inv_cum.value(unit_interval[i]);
      y[i] = f.value(x[i]);
    }
    return function<I>{x,y};
  }
  
  double significant_digits(double x, size_t n) {
    std::stringstream ss;
    ss << std::setprecision(n) << x;
    return std::stod(ss.str());
  }
  
  template<class I>
  function<I> significant_digits(function<I>& f, size_t nx, size_t ny) {
    stdvec x = f.x();
    stdvec y = f.y();
    function<I> g;
    g.header(f.header());
    for (size_t i = 0; i < x.size(); ++i) {
      double xi = significant_digits(x[i],nx);
      double yi = significant_digits(y[i],ny);
      g.append({xi,yi});
    }
    return g;
  }
  
  using pl_function = function<piecewise_linear>;
  using pp_function = function<piecewise_power>;
  using pe_function = function<piecewise_exponential>;
}

#endif
