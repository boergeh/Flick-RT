#ifndef flick_receiver
#define flick_receiver

#include "radiation_package.hpp"
#include "../numeric/histogram.hpp"

namespace flick {
  class receiver
  {
    std::vector<radiation_package> rps_;
    bool is_active_{false};
  public:
    void clear() {
      rps_.clear();
      rps_.shrink_to_fit();
    }
    void receive(const radiation_package& rp) {
      if (is_active_)
	rps_.emplace_back(rp);
    }
    size_t received_packages() {
      return rps_.size();
    }
    void activate() {
      is_active_ = true;
    }
    double radiant_flux() {
      double rf = 0;
      for(size_t i = 0; i < rps_.size(); ++i)
	rf += rps_[i].stokes().I();
      return rf;
    }
    double radiance(const unit_vector& direction, double acceptance_angle) {
      double sum = 0;
      for(size_t i = 0; i < rps_.size(); ++i) {
	unit_vector prop_dir = rps_[i].pose().z_direction(); 
	double mu = dot(prop_dir,direction);
	if (mu > cos(acceptance_angle)) {
	  sum += rps_[i].stokes().I() / fabs(prop_dir.mu());
	}
      }
      double solid_angle = 2*constants::pi*(1-cos(acceptance_angle));
      double L = sum /solid_angle;
      return L;
    }
    double mean_traveling_length() {
      double l = 0;
      for(size_t i=0; i<rps_.size(); ++i)
	l += rps_[i].traveling_length();
      return l/rps_.size();
    }
    histogram polar_angle_distribution(size_t n_emitter, size_t n_receiver) {
      double epsilon = 1e-9;
      equal_bins x_bins{-epsilon, 1+epsilon, n_emitter};
      equal_bins y_bins{-epsilon, 1+epsilon, n_receiver};
      histogram h(x_bins,y_bins);
      for(size_t i = 0; i < rps_.size(); ++i) {
	double x = rps_[i].emission_direction().mu();
	unit_vector surf_normal = {0,0,1};
	double y = dot(rps_[i].pose().z_direction(),surf_normal);
	h.add(x,y,rps_[i].stokes().I());
      }
      return h;
    }
    friend std::ostream& operator<<(std::ostream &os, const receiver& re) {
      for (auto& rp : re.rps_) {
	os << rp << '\n';
      }
      return os;
    }    
  };
}

#endif
