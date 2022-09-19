#ifndef flick_ordinary_mc
#define flick_ordinary_mc

#include "wall_interactor.hpp"
#include "material_interactor.hpp"
#include "../material/material.hpp"

namespace flick {
  using navigator = geometry::navigator<content>;
namespace transporter {
  class ordinary_mc {
    geometry::volume<content> outer_volume_;
    uniform_random rnd_;
    geometry::navigator<content> nav_;
    radiation_package rp_;
    //double sampling_asymmetry_factor_{0.8};
  public:
    ordinary_mc(const geometry::volume<flick::content>& outer_volume)
      : outer_volume_{outer_volume} {
      nav_ = geometry::navigator<flick::content>(outer_volume_);
    }
    //void set(const sampling_asymmetry_factor& saf) {
    //  sampling_asymmetry_factor_ = saf();
    //}
    bool lost_in_space() {
      return (!nav_.current_volume().has_outer_volume()
	      && !nav_.next_intersection(rp_.pose()).has_value());
    }
    flick::content& content(const std::string& volume_name) {
      return nav_.find(volume_name).content();
    }
    receiver& inward_receiver(const std::string& volume_name) {
      return nav_.find(volume_name).content().inward_receiver();
    }
    receiver& outward_receiver(const std::string& volume_name) {
      return nav_.find(volume_name).content().outward_receiver();
    }
    void transport_radiation(emitter em,
			     const std::string& emitter_volume_name,
			     double sampling_asymmetry_factor = 0.8) {
      geometry::volume<flick::content>* ev = &nav_.find(emitter_volume_name);
      while (!em.is_empty()) {
	nav_.go_to(*ev);
	rp_ = em.emit();
	double scattering_optical_depth = -log(rnd_(0,1));
	while (!rp_.is_empty() && !lost_in_space()) {
	  if (!nav_.current_volume().content().has_material()) {
	    nav_.current_volume().content().fill<material::vacuum>();
	  }
	  material::base& material = nav_.current_volume().content().material();
	  material_interactor mi(rp_,material,rnd_,scattering_optical_depth,
				 sampling_asymmetry_factor);
	  double dw = distance_to_wall();
	  double ds = mi.distance_to_scattering();
	  if (ds < dw) {
	    mi.deposite_energy_to_heat(ds);
	    mi.scatter();
	    scattering_optical_depth = -log(rnd_(0,1));
	  } else {
	    wall_interactor wi(nav_,rp_,rnd_);
	    mi.deposite_energy_to_heat(dw);
	    wi.interact_with_wall(); 
	    scattering_optical_depth -= material.scattering_optical_depth(dw);
	  }
	}
      }
    }
  private:
    double distance_to_wall() {
      std::optional<pose> p = nav_.next_intersection(rp_.pose());
      if (!p.has_value())
	return 0;
      return norm((*p).position()-rp_.pose().position()); 
    }
  };
}
}

#endif
