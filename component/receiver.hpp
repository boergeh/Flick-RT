#ifndef flick_receiver
#define flick_receiver

#include "radiation_package.hpp"
//#include "histogram.hpp"

namespace flick {
  class receiver
  {
    std::vector<radiation_package> rps_;
    //std::vector<histogram> stokes(4);
    //histogram traveling_length_;
    //histogram scattering_events_;
    bool is_active_{false};
    //bool is_inward_accepting_{true};
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
      for(size_t i=0; i<rps_.size(); ++i)
	rf += rps_[i].stokes().I();
      return rf;
    }
    double radiant_intensity(const unit_vector& direction,
			     double acceptance_angle) {
      double mu_max = 2*constants::pi*(1-cos(acceptance_angle));
      double sum_ri = 0;
      double sum_omega = 0;
      for(size_t i = 0; i < rps_.size(); ++i) {
	double mu = dot(rps_[i].pose().z_direction(),direction);
	if (mu < mu_max) {
	  sum_ri += rps_[i].stokes().I(); 
	  sum_omega += 2*constants::pi*(1-mu);
	}
      }
      return sum_ri/sum_omega;
    }
    double mean_traveling_length() {
      double l = 0;
      for(size_t i=0; i<rps_.size(); ++i)
	l += rps_[i].traveling_length();
      return l/rps_.size();
    }
    //bool is_active() {
    //  return is_active_;
    //}
    //void outward_accepting() {
    //  is_inward_accepting_ = false;
    //}
    // if size larger than so and so fill histogram empty rps.
    /*
    histogram stokes(enum element, size_t n_bins) const {
      histrogram h(equal_bins{-pi,pi,n_bins})
      for (size_t i=0; i < rps_.size(); ++i) {
	h.add()
      }
      
    }
    */
    friend std::ostream& operator<<(std::ostream &os, const receiver& re) {
      for (auto& rp : re.rps_) {
	os << rp << '\n';
      }
      return os;
    }    
  };
}

#endif