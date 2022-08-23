#ifndef flick_ordinary_mc
#define flick_ordinary_mc

#include "wall_interactor.hpp"

namespace flick {
namespace transporter {
  class material_interactor {
    radiaiton_package& rp_;
    material& m_;
    uniform_random& ur_;
    double scattering_optical_depth_;
    double g_;
    unit_vector scattering_direction_;
  public:
    material_interactor(radiation_package& rp,
			material& m,
			const uniform_random& ur,
			double scattering_optical_depth,
			double sampling_asymmetry_factor)
      : rp_{rp}, m_{m}, ur_{ur}, g_{sampling_asymmetry_factor},
	scattering_optical_depth_{scattering_optcal_depth} {
      m_.position(rp.position());
      m_.direction(rp.direction());
      distance_to_scattering_ = m_.column_length(scattering_optical_depth_);
    }
    void move_to_scattering_event() {
      rp.move(distance_to_scattering_);
    }
    void deposite_energy_to_heat() {
      rp_.scale_intensity(exp(-m_.optical_depth(distance_to_scattering_)));
    }
    double distance_to_scattering() {
      return distance_to_scattering_;
    }
    void find_scattering_direction() {
      double& g = sampling_asymmetry_factor_;
      double polar_angle = henyey_greenstein{g}.inverted_accumulated_angle(ur_(0,1));
      double azimuth_angle = ur_(0,2*constants::pi);
      pose p = rp_.pose();
      p.rotate_about_z(azimuth_angle);
      p.rotate_about_x(polar_angle);
      scattering_direction_ =  p.z_direction();
    }
    void reorient_traveling_direction() {
      rp_.traveling_direction(scattering_direction_);
    }
    void reshape_polarization() {
      set_x_axis_parallel_with_scattering_plane();
      m_.direction(scattering_direction_);
      rp_.reshape_polarization(m_.mueller_matrix());
    }
    void likelihood_scale_intensity() {
      double hg = henyey_greenstein{g}.value(scattering_angle());
      m_.direction(new_direction);
      rp_.scale_intensity(m_.mueller_matrix().value(0,0)/hg);
    }
  private:
    double scattering_angle() {
      return acos(dot(rp_.pose().z_direction(),scattering_direction_));
    }
    void set_x_axis_parallel_with_scattering_plane() {
      pose& p = rp_.pose();
      vector n = cross(p.z_direction(),scattering_direction_);
      double ang = acos(dot(n, p.x_direction()));
      rp_.rotate_about_local_z(ang);
    }

  };
  class ordinary_mc {
    geometry::volume<content>& outer_volume_;
    geometry::volume<content>& emitter_volume_;
    emitter& em_;
    uniform_random ur_;
    coating::lambert_angle_generator lag_;
  public:    
    ordinary_mc(geometry::volume<content>& outer_volume,
		geometry::volume<content>& emitter_volume,
		emitter& em)
      : outer_volume_{outer_volume}, emitter_volume_{emitter_volume},
	em_{em} {}
    bool lost_in_space(geometry::navigator<content>& nav,
		       radiation_package& rp) {
      if (!nav.current_volume().has_outer_volume()
	  && rp.weighted_traveling_length() > 0)
	return true;
      return false;
    }
    void move_to_event() {
      
    }
    void run() {
      auto nav =  geometry::navigator<content>(outer_volume_);
      while (!em_.is_empty()) {
	nav.go_to(emitter_volume_);
	radiation_package rp = em_.emit();
	bool has_scattered = true;
	while (!rp.is_empty() && !lost_in_space(nav,rp)) {
	  if (has_scattered)
	    scattering_optical_depth = -log(ur_()); //include likelihood?
	  material::material& m = nav.current_volume().content().material();
	  material_interactor mi(rp,m,ur_,scattering_optical_depth,0.8);
	  wall_interactor wi(nav,rp,ur_,lag_);
	  if (mi.distance_to_scattering() < wi.distance_to_wall()) {
	    mi.move_to_scattering_event();
	    mi.deposite_energy_to_heat();
	    mi.find_scattering_direction();
	    mi.likelihood_scale_intensity();
	    mi.reorient_traveling_direction();
	    mi.reshape_polarization();
	    has_scattered = true;
	  } else if (wi.is_reflected()) {
	    consume_scattering_optical_depth(ma);
	    wi.move_close_to_intersection();
	    //wi.move_small_step_backward();
	    wi.deposite_energy_to_heat();
	    //wi.find_new_orientation();
	    wi.reorient_traveling_direction();
	    wi.reshape_polarization();
	    wi.likelihood_scale_intensity();
	  } else if (wi.is_transmitted()) {
	    wi.move_through_wall();
	    //wi.small_step_forward();
	    wi.increase_activated_receiver();
	    wi.reorient_traveling_direction();
	    wi.reshape_polarization();
	    wi.likelihood_scale_intensity();
	  } else if (wi.is_absorbed_by_wall()) {
	    wi.extinguish_radiation();
	  }  
	}
      }
    }
  };
}
}

#endif
