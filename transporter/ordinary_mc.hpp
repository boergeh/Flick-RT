#ifndef flick_ordinary_mc
#define flick_ordinary_mc

#include "mc_basic.hpp"

namespace flick {
namespace transporter {
  class ordinary_mc {
    //geometry::volume<content>* v;
  public:
    number_of_packages np_;
    ordinary_mc(const number_of_packages& np)
      : np_{np} {}
    //      basic_mc(const precision& p) {	
    //}
    //tbi
    
    //void run(geometry::volume<content>& v,
    //       geometry::volume<content>& emitter_volume,
    //       emitter& em) {
    //mc_basic(v).run(em, emitter_volume);
  };
}
}

#endif
