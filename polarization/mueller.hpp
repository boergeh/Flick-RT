#ifndef flick_mueller
#define flick_mueller
#include "../numeric/function.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/std_operators.hpp"


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

  std::vector<double> hg_importance_sampling(double g, size_t n_points) {
    henyey_greenstein hg(g);
    std::vector<double> fraction = range(0,1,n_points).linspace();
    std::vector<double> x(n_points); 
    for (size_t i=0; i<n_points; ++i) {
      x[i] = hg.inverted_accumulated_angle(fraction[i]); 
    }
    return x;
  }

  template<class Function>
  Function resample(const Function& f, const std::vector<double>& new_x) {
    std::vector<double> new_y(new_x.size());
    for (size_t i=0; i<new_x.size(); ++i) {
      new_y[i] = f.value(new_x[i]);
    }
    return Function{new_x, new_y};
  }

  template<class Function>
  Function to_cos_x(const Function& f) {
    Function fc;
    size_t n = f.size();
    for (size_t i=0; i<n; ++i) {
      fc.append({cos(f.x()[n-i-1]), f.y()[n-i-1]});
    }
    return fc;
  }
  
  template<class Function>
  Function to_acos_x(const Function& f) {
    Function fc;
    size_t n = f.size();
    for (size_t i=0; i<n; ++i) {
      fc.append({acos(std::clamp<double>(f.x()[n-i-1],-1,1)), f.y()[n-i-1]});
    }
    return fc;
  }
  
  template<class Function>
  Function hg_phase_function(double g, size_t n_points,
			     double sampling_skew_factor = 0.8) {
    henyey_greenstein hg(g);
    std::vector<double> x = hg_importance_sampling(sampling_skew_factor, n_points);
    x = vec::cos(x);
    std::vector<double> y(x.size());
    for (size_t i=0; i<x.size(); ++i) {
      y[i] = hg.value(x[i]);
    }
    return Function{x,y};    
  }
 
  class tabulated_phase_function {
    pe_function p_;
    pl_function residual_;
  public:
    tabulated_phase_function() = default;
    tabulated_phase_function(const pe_function& p)
      : p_{p} {
      double epsilon = 1e-6;
      ensure(p.x().front() > -1-epsilon && p.x().back() < 1+epsilon);
    }  
    const std::vector<double>& x() const {
      return p_.x();
    }    
    const std::vector<double>& y() const {
      return p_.y();
    }    
    bool empty() const {
      if (p_.size()==0)
	return true;
      return false;
    }
    double value(double mu) const {
      return p_.value(mu);
    }
    double asymmetry_factor() const
    // Represent phase function as Henyey-Greenstein plus a residual for
    // more accurate interpolation and integration.
    {
      double g = 0.9;
      henyey_greenstein hg(g);
      pl_function residual;
      for (size_t i=0; i<p_.size();i++)
	residual.append({p_.x()[i], p_.y()[i] - hg.value(p_.x()[i])});
      return g + 2*constants::pi*pl_function{residual.x(),
	residual.x()*residual.y()}.integral(-1,1);
    }
    double integral_4pi() const {
      return 2*constants::pi*p_.integral(-1,1);
    }  
  private:
    void ensure(bool b) const {
      if (!b)
	throw std::runtime_error("polarization mueller tabulated_phase_function");
    }
  };
   
  class angular_mueller
  // Keeps Mueller matrix values at discrete angles 
  {
    tabulated_phase_function p_;
    struct element {
      size_t row;
      size_t col;
      pl_function f;
    };
    std::vector<element> elements_;
    double epsilon = 1e-6;
  public:
    angular_mueller() = default;
    angular_mueller(const tabulated_phase_function& p)
      : p_{p} {
      //pe_function p2 = to_acos_x(pe_function{p.x(),p.y()});
      //std::cout << p2;
      // add(0,0,pl_function{p2.x(),p2.y()});
      add(0,0,pl_function{p.x(),p.y()});
    }
    angular_mueller& append(size_t row, size_t col,
			    double mu, double value) {
      ensure(mu > -1-epsilon && mu < 1+epsilon);
      auto& e = elements_;
      for (size_t i=0; i<e.size(); ++i) {
	if (e[i].row == row && e[i].col == col) {
	  e[i].f.append({mu,value});
	  return *this;
	}
      }
      element new_e;
      new_e.row = row;
      new_e.col = col;
      new_e.f = pl_function{};
      e.push_back(new_e);
      e.back().f.append({mu,value});
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
    double value(size_t row, size_t col, double mu) const {
      ensure(mu > -1-epsilon && mu < 1+epsilon);
      auto& e = elements_;
      for (size_t n = 0; n < e.size(); ++n) {
	if (e[n].row == row && e[n].col == col) {
	  return e[n].f.value(mu)*p_.value(mu);
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
      update_phase_function(am, w);
      return *this;
    }
    void print(double mu) {
      auto& e = elements_;
      std::cout << std::endl;
      for (size_t i = 0; i < 4; ++i) {
	for (size_t j = 0; j < 4; ++j) {
	  std::cout << value(i,j,mu) << " ";
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
      pe_function p = to_acos_x(pe_function{p_.x(),p_.y()});
      std::vector<double> x = p.x();
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
      pe_function p = pe_function{p_.x(),p_.y()};
      pe_function amp = pe_function{am.p_.x(),am.p_.y()};
      pl_function p1{p.x(),p.y()};
      pl_function p2{amp.x(),amp.y()};
      std::vector<double> x = p1.x();
      if (p2.size() > p1.size())
	x = p2.x();
      std::vector<double> y(x.size());
      for (size_t i=0; i<x.size(); ++i) {
	y[i] = (1-w)*p1.value(x[i]) + w*p2.value(x[i]);
      }      
      //p_ = tabulated_phase_function(to_cos_x(pe_function{x,y}));
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
    void ensure(bool b) const {
      if (!b)
	throw std::runtime_error("polarization angular_mueller");
    }
  };
}

#endif
