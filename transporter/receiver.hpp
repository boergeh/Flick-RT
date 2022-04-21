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
  public:
    void receive(const radiation_package& rp) {
      rps_.emplace_back(rp);
    }
    size_t received_packages() {
      return rps_.size();
    }
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
