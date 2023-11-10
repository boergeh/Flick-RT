#ifndef flick_mueller
#define flick_mueller
#include "../numeric/function.hpp"
#include "../numeric/physics_function.hpp"

namespace flick {
  class mueller
  // Holds non-zero Mueller matrix element values for one angle
  {
    struct element {
      size_t row;
      size_t col;
      double value;
    };
    std::vector<element> elements_;
  public:
     mueller& add(size_t row, size_t col, double value) {
      elements_.push_back(element{row,col,value});
      return *this;
    }
    size_t size() const {
      return elements_.size();
    }
    const element& operator()(size_t n) const {
      return elements_.at(n);
    }
    double value(size_t row, size_t col) const {
      auto& e = elements_;
      for (size_t n = 0; n < e.size(); ++n) {
	if (e[n].row == row && e[n].col == col) {
	  return e[n].value;
	}
      }
      return 0;
    }
  };

  std::ostream& operator<<(std::ostream &os, const mueller& m) {
    size_t n = 0;
    os << '\n';
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
	if (n < m.size() && m(n).row == i && m(n).col == j) {
	  os << m(n).value << " ";
	  n++;
	} else {
	  os << 0 << " ";
	}
      }
      os << '\n';
    }
    return os;
  }

  class tabulated_phase_function
  // Represent tabulated phase function as a Henyey-Greenstein with
  // residuals.
  {
    double g_;
    pl_function residuals_;
    pe_function p_mu_;
  public:
    tabulated_phase_function() = default;
    tabulated_phase_function(const pe_function& p) {
      if (p.x()[0] < 0)
	p_mu_ = p;
      else
	cosine_transform_angles(p);	 
      g_ = asymmetry_factor_(p_mu_);
      for (size_t i=0; i<4; ++i) {
	set_residuals(p_mu_);
	pe_function pe_residuals{residuals_.x(),residuals_.y()};
	g_ += asymmetry_factor_(pe_residuals);
      }
      set_residuals(p_mu_);
    }
    std::vector<double> x() const {
      std::vector<double> new_x(residuals_.size());
      for (size_t i=0; i<new_x.size(); ++i)
	new_x[new_x.size()-i-1] = acos(residuals_.x()[i]);
      return new_x;
    }
    std::vector<double> y() const {
      std::vector<double> xx = x();
      std::vector<double> new_y(xx.size());
      for (size_t i=0; i<xx.size(); ++i)
	new_y[i] = value(xx[i]);
      return new_y;
    }
    bool empty() {
      if (residuals_.size()==0)
	return true;
      return false;
    }
    double value(double angle) const {
      return flick::henyey_greenstein(g_).phase_function(angle)
	+ residuals_.value(cos(angle));
    }
    double asymmetry_factor() const {
      return g_;
    }
    double integral() const {
      return 1/(2*constants::pi) + residuals_.integral(-1,1);
    }
  private:
    void cosine_transform_angles(const pe_function& p) {
      size_t n = p.size();
      for (size_t i=0; i<n; ++i) {
	double mu = cos(p.x()[n-i-1]);
	p_mu_.append({mu,p.y()[n-i-1]});
      }
    }
    void set_residuals(const pe_function& p) {
      residuals_.clear();
      for (size_t i=0; i<p.size(); ++i) {
	double mu = p.x()[i];
	double r = p.y()[i]-flick::henyey_greenstein(g_).
	  phase_function(acos(mu));
	residuals_.append({mu, r});
      }
    }
    double asymmetry_factor_(const pe_function& p) {
      pl_function p2;
      for (size_t i=0; i<p.size(); ++i) {
	double mu = p.x()[i];
	p2.append({mu, p.y()[i]*mu});
      }
      return 2*constants::pi*p2.integral(-1,1);
    }    
  };

  pe_function hg_phase_function(double g, size_t n_points) {
    auto x = range(0,constants::pi,n_points).linspace();
    std::vector<double> y(x.size());
    for (size_t i=0; i<y.size(); ++i)
      y[i] = henyey_greenstein(g).phase_function(x[i]);
    return pe_function{x,y};    
  }
 
  class angular_mueller {
    tabulated_phase_function p_;
    struct element {
      size_t row;
      size_t col;
      pl_function f;
    };
    std::vector<element> elements_;
  public:
    angular_mueller() = default;
    angular_mueller(const tabulated_phase_function& p)
      : p_{p} {
      add(0,0,pl_function{p.x(),p.y()});
    }
    angular_mueller& append(size_t row, size_t col,
			    double angle, double value) {
      auto& e = elements_;
      for (size_t i=0; i<e.size(); ++i) {
	if (e[i].row == row && e[i].col == col) {
	  e[i].f.append({angle,value});
	  return *this;
	}
      }
      element new_e;
      new_e.row = row;
      new_e.col = col;
      new_e.f = pl_function{};
      e.push_back(new_e);
      e.back().f.append({angle,value});
      return *this;
    }
    angular_mueller& add(size_t row, size_t col,
			 const pl_function& f) {
      elements_.push_back(element{row,col, normalize(f)});
      return *this;
    }
    size_t size() const {
      return elements_.size();
    }
    const element& operator()(size_t n) const {
      return elements_.at(n);
    }
    double value(size_t row, size_t col, double angle) const {
      auto& e = elements_;
      for (size_t n = 0; n < e.size(); ++n) {
	if (e[n].row == row && e[n].col == col) {
	  return e[n].f.value(angle)*p_.value(angle);
	}
      }
      return 0;
    }
    angular_mueller& add(const angular_mueller& am, double fractional_weight) {
      if (p_.empty()) {
	*this = am; 
	return *this;
      }
      double w = fractional_weight;
      auto& e1 = elements_;
      auto& e2 = am.elements_;
      std::vector<element> new_e;
      std::vector<bool> added1(e1.size(),false);
      std::vector<bool> added2(e2.size(),false);
      for (size_t i=0; i<e1.size(); ++i) {
	for (size_t j=0; j<e2.size(); ++j) {
	  if (e1[i].row == e2[j].row && e1[i].col == e2[j].col) {
	    added1[i]=true;
	    added2[j]=true;
	    element e;
	    e.row = e1[i].row;
	    e.col = e1[i].col;
	    e.f = add_elements(e1[i].f, p_, e2[j].f, am.p_, w);
	    new_e.push_back(e);
	  }
	}
      }
      for (size_t i=0; i<e1.size(); ++i) {
	if (!added1[i]) {
	  element e;
	  e.row = e1[i].row;
	  e.col = e1[i].col;
	  e.f = add_elements(e1[i].f, p_, {0}, am.p_, w);
	  new_e.push_back(e);
	}
      }
      for (size_t i=0; i<e2.size(); ++i) {
	if (!added2[i]) {
	  element e;
	  e.row = e2[i].row;
	  e.col = e2[i].col;
	  e.f = add_elements({0}, p_, e2[i].f, am.p_, w);
	  new_e.push_back(e);
	}	  
      }
      elements_ = new_e;
      update_phase_function(am,fractional_weight);
      return *this;
    }
    void print(double angle) {
      auto& e = elements_;
      std::cout << std::endl;
      for (size_t i = 0; i < 4; ++i) {
	for (size_t j = 0; j < 4; ++j) {
	  std::cout << value(i,j,angle) << " ";
	}
	std::cout << std::endl;
      }
    }
  private:
    pl_function add_elements(const pl_function& f1,
			     const tabulated_phase_function& p1,
			     const pl_function& f2,
			     const tabulated_phase_function& p2,
			     double w) {
      std::vector<double> x = p_.x();
      std::vector<double> y(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	double p = (1-w)*p1.value(x[i]) + w*p2.value(x[i]);
	y[i] = ((1-w)*f1.value(x[i])*p1.value(x[i]) +
		w*f2.value(x[i])*p2.value(x[i])) / p;
      }
      return pl_function{x,y};
    }
    void update_phase_function(const angular_mueller& am,
				double fractional_weight) {
      double w = fractional_weight;
      pl_function p1{p_.x(),p_.y()}; 
      pl_function p2{am.p_.x(),am.p_.y()};
      std::vector<double> x = p1.x();
      if (p2.size()>p1.size())
	x = p2.x();
      std::vector<double> y(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	y[i] = (1-w)*p1.value(x[i]) + w*p2.value(x[i]);
      }
      p_ = tabulated_phase_function(pe_function{x,y});
    }
    pl_function normalize(const pl_function& f) {
      pl_function f_new;
      for (size_t i=0; i < f.size(); ++i) {
	double x = f.x()[i];
	double y = f.y()[i] / p_.value(x);
	f_new.append({x,y});
      }
      return f_new;
    }
  };
}

#endif
