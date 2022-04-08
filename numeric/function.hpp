#ifndef flick_function
#define flick_function

#include <optional>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "sorted_vector.hpp"

namespace flick {
  struct point {
    double x_;
    double y_;
  public:
    point(){}
    point(double x, double y) : x_(x), y_(y){}
    double x() const {return x_;}
    double y() const {return y_;}
  };
 
  class basic_interpolation {
  protected:
    point p1;
    point p2;
    double a;
    double b;
  public:
    basic_interpolation(const point& low, const point& high)
      : p1{low}, p2{high} {
    }
  };
  
  class piecewise_linear : public basic_interpolation {
  public:
    static step_type step_type() {
      return step_type::linear;
    }
    piecewise_linear(const point& low, const point& high)
      : basic_interpolation(low, high) {
      a = (p2.y()-p1.y())/(p2.x()-p1.x());
      b = p1.y()-a*p1.x();
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
      if (fabs(a) < std::numeric_limits<double>::epsilon()) {
	if (fabs(b) < std::numeric_limits<double>::epsilon())
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
  public:
    static step_type step_type() {
      return step_type::linear;
    }
    piecewise_exponential(const point& low, const point& high)
      : basic_interpolation(low, high) {
      ensure(low.y() > 0 && high.y() > 0);
      b = (log(p2.y())-log(p1.y())) / (p2.x()-p1.x());
      a = exp(log(p1.y()) - b * p1.x());
    }
    double y(double x) {
      return a * exp(b*x);
    }
    double derivative(double x) {
      return a * b * exp(b*x);
    }
    double integral(double limit_a, double limit_b) {
      if (fabs(b) < std::numeric_limits<double>::epsilon())
	return a*(limit_b-limit_a);
      return a/b*(exp(b*limit_b)-exp(b*limit_a));
    }
    std::optional<double> integral_limit_b(double limit_a,
					   double integral_value) {
      if (fabs(b) < std::numeric_limits<double>::epsilon())
	return limit_a + integral_value / a;
      double arg = integral_value*b/a + exp(b*limit_a);
      if (arg <= 0)
	return std::nullopt;
      return log(arg)/b;
    }
  private:
    void ensure(bool b) {
      if (!b) {
	throw std::invalid_argument("piecewise_exponential");
      }
    }    
  };

  class piecewise_power : public basic_interpolation {
  public:
    static step_type step_type() {
      return step_type::exponential;
    }
    piecewise_power(const point& low, const point& high)
      : basic_interpolation(low, high) {
      ensure(low.x() > 0 && high.x() > 0);
      ensure(low.y() > 0 && high.y() > 0);
      b = (log(p2.y())-log(p1.y())) / (log(p2.x())-log(p1.x()));
      a = exp(log(p1.y()) - b * log(p1.x()));
    }
    double y(double x) {
      return a*pow(x,b);
    }
    double derivative(double x) {
      if (fabs(b) < std::numeric_limits<double>::epsilon())
	return 0;
      return a / b * pow(x,b-1);
    }
    double integral(double limit_a, double limit_b) {
      ensure(limit_a > 0 && limit_b > 0);
      if (fabs(b+1) < std::numeric_limits<double>::epsilon())
	return a*log(limit_b/limit_a);
      return a/(b+1)*(pow(limit_b,b+1)-pow(limit_a,b+1));
    }
    std::optional<double> integral_limit_b(double limit_a,
					   double integral_value) {
      if (fabs(b+1) < std::numeric_limits<double>::epsilon())
	return limit_a * exp(integral_value/a);
      double arg = integral_value*(b+1)/a + pow(limit_a,b+1);
      if (arg <= 0)
	return std::nullopt;
      return pow(arg,1/(b+1));
    }
  private:
    void ensure(bool b) {
      if (!b) {
	throw std::invalid_argument("piecewise_power");
      }
    }    
  };
 
  template<class Interpolation>
  class function
  {
    std::string header_;
    std::vector<double> yv_;
    sorted_vector xv_;
    std::shared_ptr<sorted_vector::iterator> it;
  public:
    function() {};
    function(double value) : xv_{{1}}, yv_{std::vector<double>{value}} {}
    function(const std::vector<double>& xv, const std::vector<double>& yv)
      : xv_{xv}, yv_{yv} {
      ensure(xv.size()==yv.size() && xv.size() > 1);
      xv_.set_step_type(Interpolation::step_type());
    }
    void add_extrapolation_points(double weight=1) {
      ensure(xv_.size() > 1);
      if (Interpolation::step_type()!=step_type::linear)
	ensure(weight > 0);
      std::vector<double> xv = xv_.all_values();
      double dx_front = xv[1]-xv[0];
      double dx_back = xv[xv.size()-1]-xv[xv.size()-2];
      if (Interpolation::step_type()!=step_type::linear) {
	dx_front = xv[0]*dx_front/xv[1];
	dx_back = xv.back()*dx_back/xv.back();
      }
      xv.insert(xv.begin(),xv[0]-dx_front/2);
      xv.insert(xv.begin(),xv[0]-dx_front);
      xv.emplace_back(xv.back()+dx_back/2);
      xv.emplace_back(xv.back()+dx_back);
      xv_ = sorted_vector{xv};
      yv_.insert(yv_.begin(),yv_[0]*weight);
      yv_.insert(yv_.begin(),yv_[0]*weight);
      yv_.emplace_back(yv_.back()*weight);
      yv_.emplace_back(yv_.back()*weight);
    }
    void header(const std::string& h) {
      header_ = h;
    }
    std::string header() const {
      return header_;
    }
    void append(const point& p) {
      xv_.append(p.x());
      yv_.push_back(p.y());
    }
    void scale_x(double factor) {
      ensure(factor > 0);
      xv_.scale(factor);
    }
    void scale_y(double factor) {
      ensure(factor > 0);
      for (size_t i=0; i<yv_.size(); ++i)
	yv_[i] *= factor;
    }
    std::vector<double> x() {
      return xv_.all_values();
    }
    std::vector<double> y() {
      return yv_;
    }
    double value(double x=1) {
      if (yv_.size()==1)
	return yv_[0];
      ensure(yv_.size() > 1);
      auto [p1, p2] = points_at(x);
      return Interpolation{p1,p2}.y(x);      
    }
    double derivative(double x) {
      ensure(yv_.size() > 1);
      auto [p1, p2] = points_at(x);
      return Interpolation{p1,p2}.derivative(x);
    }
    std::optional<double> integral_limit_b(double limit_a,
					   double integral_value) {
      ensure(yv_.size() > 1);
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
	double next_area = Interpolation{p1,p2}.integral(limit_a, p2.x());
	if(fabs(area + next_area) > fabs(integral_value) || it->is_in_end_bin())
	  return Interpolation{p1,p2}.integral_limit_b(limit_a,
						       integral_value-area);
	area += next_area;
	limit_a = p2.x();
	it->move_to_next_bin();
      }      
    }
    double integral(double limit_a, double limit_b) {      
      ensure(yv_.size() > 1);
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
	if(fabs(p2.x()) > fabs(limit_b) || it->is_in_end_bin())
	  return area + Interpolation{p1,p2}.integral(limit_a, limit_b);
	area += Interpolation{p1,p2}.integral(limit_a, p2.x());
	limit_a = p2.x();
	it->move_to_next_bin();
      }      
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const function<Interpolation>& f) {
      os << f.header();
      for (size_t i = 0; i<f.yv_.size(); ++i)
	os << f.xv_[i] << " " << f.yv_[i] << '\n';
      return os;
    }
    friend std::istream& operator>>(std::istream &is,
				    function<Interpolation>& f) {
      f.header(f.read_header(is));
      double x, y;
      while(is >> x >> y)
	f.append({x,y});
      return is;
    }
  private:
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
    point next_point(sorted_vector::iterator *it) {
      return point{xv_[it->next_index()],yv_[it->next_index()]};
    }
    point previous_point(sorted_vector::iterator *it) {
      return point{xv_[it->previous_index()],yv_[it->previous_index()]};
    }
    std::tuple<point, point> points_at(double x) {
      size_t n = xv_.find(x);
      point p1{xv_[n],yv_[n]};
      point p2{xv_[n+1],yv_[n+1]};
      return {p1,p2};
    }
    void ensure(bool b) {
      if (!b)
	throw std::invalid_argument("function");
    }
  };
  /*
  template<class Interpolation>
  function<Interpolation> add_extrapolation_points(function<Interpolation>& f,
						   double dx_factor=0.05) {
    std::vector<double> x = f.x();
    std::vector<double> y = f.y();
    function<Interpolation> g;
    g.header(f.header());
    g.append({x.at(0)*(1-dx_factor),y.at(0)});
    for (size_t i = 0; i < x.size(); ++i) {
      g.append({x[i],y[i]});
    }
    g.append({x.back()*(1+dx_factor),y.back()});
    return g;
  }
  */
  double significant_digits(double x, size_t n) {
    std::stringstream ss;
    ss << std::setprecision(n) << x;
    return std::stod(ss.str());
  }
  template<class Interpolation>
  function<Interpolation> significant_digits(function<Interpolation>& f,
					     size_t nx, size_t ny) {
    std::vector<double> x = f.x();
    std::vector<double> y = f.y();
    function<Interpolation> g;
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
