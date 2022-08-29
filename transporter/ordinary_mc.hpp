#ifndef flick_ordinary_mc
#define flick_ordinary_mc

#include "wall_interactor.hpp"
#include "material_interactor.hpp"
#include "../material/material.hpp"

namespace flick {
namespace transporter {  
  class ordinary_mc {
    geometry::volume<content>& outer_volume_;
    geometry::volume<content>& emitter_volume_;
    emitter& emitter_;
    uniform_random random_;
    coating::lambert_angle_generator lag_;
    geometry::navigator<content> nav_;
    radiation_package rp_;
  public:    
    ordinary_mc(geometry::volume<content>& outer_volume,
		geometry::volume<content>& emitter_volume,
		emitter& em)
      : outer_volume_{outer_volume}, emitter_volume_{emitter_volume},
	emitter_{em} {
      nav_ = geometry::navigator<content>(outer_volume_);
    }
    bool lost_in_space() {
      return (!nav_.current_volume().has_outer_volume()
	      && !nav_.next_intersection(rp_.pose()).has_value());
    }
    void run(double sampling_optical_depth=1,
	     double sampling_asymmetry_factor=0.8) {     
      material::vacuum vacuum;
      while (!emitter_.is_empty()) {
	nav_.go_to(emitter_volume_);
	rp_ = emitter_.emit();
	double scattering_optical_depth = -log(random_(0,1));
	while (!rp_.is_empty() && !lost_in_space()) {
	  if (!nav_.current_volume().content().has_material()) {	  
	    nav_.current_volume().content().fill<material::vacuum>();
	  }
	  material::base& material = nav_.current_volume().content().material();
	  material_interactor mi(rp_,material,random_,scattering_optical_depth,
				 sampling_asymmetry_factor);
	  wall_interactor wi(nav_,rp_,random_,lag_);
	  if (mi.distance_to_scattering() < wi.distance_to_wall()) {
	    mi.deposite_energy_to_heat(mi.distance_to_scattering());
	    mi.scatter();
	    scattering_optical_depth = -log(random_(0,1));
	    //include likelihood?
	  } else {
	    mi.deposite_energy_to_heat(wi.distance_to_wall());
	    wi.interact_with_wall();
	    scattering_optical_depth -= material.scattering_optical_depth(wi.distance_to_wall());
	  }
	  //std::cout << rp_ << std::endl;
	}
      }
    }
  };
}
}

#endif
