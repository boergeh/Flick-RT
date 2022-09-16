#ifndef flick_material_tabulated
#define flick_material_tabulated

#include "monocrome_iop.hpp"
#include "../numeric/function.hpp"

namespace flick {
namespace material {
  class phase_function {
    pl_function p_;
    pl_function residuals_;
    double g_;
    const double pi = constants::pi;
    //pp_function inverted_accumulated_;
  public:
    phase_function(const pl_function& p) : p_{p} {
      // check
      std::vector<double> x = p.x();
      std::vector<double> y = p.y();
      std::vector<double> mu(x.size());
      std::vector<double> y2(x.size());
      for (size_t i=0; i<mu.size(); ++i) {
	mu[i] = cos(x[mu.size()-i-1]);
	y2[i] = y[mu.size()-i-1];
      }
      p_ = pl_function(mu,y2);
      double g = compute_asymmetry_factor(p_);
      compute_residuals(g);
      g_ = g + compute_asymmetry_factor(residuals_);

      g = asymmetry_factor();
      compute_residuals(g);
      g_ = g + compute_asymmetry_factor(residuals_);

      //compute_residuals(g);
      //renormalize();
    }
    double value(double angle) const {
      return p_.value(angle);
    }
    double asymmetry_factor() const {
      return g_;
    }
    double integral() {
      return 1/(2*pi) + residuals_.integral();
    }
  private:
    //  void renormalize() {
    //    double k = (1-1/(2*pi))/integral();
    //    residuals_.scale_y(k);
    // }
    void compute_residuals(double g) {
      auto y_hg = compute_hg(g);
      std::vector<double> x = p_.x();
      std::vector<double> y = p_.y();
      std::vector<double> y2(x.size());
      for (size_t i=0; i<x.size(); ++i)
	y2[i] = y[i] - y_hg[i];
      residuals_ = pl_function(x,y2);
    }
    std::vector<double> compute_hg(double g) {
      std::vector<double> x = p_.x();
      std::vector<double> y(x.size());
      for (size_t i=0; i<x.size(); ++i)
	y[i] = flick::henyey_greenstein(g).phase_function(acos(x[i]));
      return y;
    }
    double compute_asymmetry_factor(const pl_function& p) {
      std::vector<double> x = p.x();
      std::vector<double> y = p.y();
      for (size_t i=0; i<y.size(); ++i) {
	y[i] *= x[i];
      }
      pl_function p2{x,y};
      return p2.integral()*2*constants::pi;
    }    
  };
  class tabulated : public monocrome_iop {
    phase_function p_;
  public:
    tabulated(const flick::absorption_coefficient& ac,
	      const flick::scattering_coefficient& sc,
	      const phase_function& p,
	      double real_refractive_index = 1)
      : monocrome_iop{ac,sc,flick::asymmetry_factor{p.asymmetry_factor()},
      real_refractive_index}, p_{p} {
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) {
      mueller m;
      double theta = angle(scattering_direction);
      m.add(0,0,p_.value(theta));
      return m;
    }
  };
}
}

#endif
