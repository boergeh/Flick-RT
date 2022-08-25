#ifndef flick_ordinary_mc
#define flick_ordinary_mc

#include "wall_interactor.hpp"
#include "../material/material.hpp"

namespace flick {
namespace transporter {
  class material_interactor {
    radiation_package& rp_;
    material::base& m_;
    uniform_random& ur_;
    double scattering_optical_depth_;
    double g_;
    unit_vector scattering_direction_;
    double scattering_angle_;
    double distance_to_scattering_;
  public:
    material_interactor(radiation_package& rp,
			material::base& m,
			uniform_random& ur,
			double scattering_optical_depth,
			double sampling_asymmetry_factor)
      : rp_{rp}, m_{m}, ur_{ur}, g_{sampling_asymmetry_factor},
	scattering_optical_depth_{scattering_optical_depth} {
      m_.set(rp_.pose());
      distance_to_scattering_ = m_.scattering_distance(scattering_optical_depth_);
      scattering_angle_ = m_.angle(scattering_direction_);
    }
    void move_to_scattering_event() {
      rp_.move(distance_to_scattering_);
    }
    void deposite_energy_to_heat(double distance) {
      rp_.scale_intensity(exp(-m_.absorption_optical_depth(distance)));
    }
    double distance_to_scattering() {
      return distance_to_scattering_;
    }
    void find_scattering_direction() {
      double polar_angle = henyey_greenstein{g_}.inverted_accumulated_angle(ur_(0,1));
      double azimuth_angle = ur_(0,2*constants::pi);
      pose p = rp_.pose();
      p.rotate_about_local_z(azimuth_angle);
      p.rotate_about_local_x(polar_angle);
      scattering_direction_ =  p.z_direction();
    }
    void reorient_traveling_direction() {
      rp_.rotate_about_local_y(scattering_angle_);
    }
    void reshape_polarization() {
      set_x_axis_parallel_with_scattering_plane();
      rp_.reshape_polarization(m_.mueller_matrix(scattering_direction_));
    }
    void likelihood_scale_intensity() {
      double hg = henyey_greenstein{g_}.phase_function(scattering_angle_);
      rp_.scale_intensity(m_.mueller_matrix(scattering_direction_).value(0,0)/hg);
    }
    void scatter() {
      move_to_scattering_event();
      find_scattering_direction();
      likelihood_scale_intensity();
      reshape_polarization();
      reorient_traveling_direction();
    }
  private:
    void set_x_axis_parallel_with_scattering_plane() {
      const pose& p = rp_.pose();
      vector n = cross(p.z_direction(),scattering_direction_);
      double ang = acos(dot(n, p.x_direction()));
      rp_.rotate_about_local_z(ang);
    }
  };
  
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
	      && rp_.weighted_traveling_length() > 0);
    }
    void run(double sampling_optical_depth, double sampling_asymmetry_factor) {     
      while (!emitter_.is_empty()) {
	nav_.go_to(emitter_volume_);
	rp_ = emitter_.emit();
	double scattering_optical_depth = -log(random_(0,1));
	while (!rp_.is_empty() && !lost_in_space()) {
	  if (!nav_.current_volume().content().has_material())	  
	    nav_.current_volume().content().fill<material::vacuum>();
	  material::base& material = nav_.current_volume().content().material();
	  //std::cout << "material set" << std::endl;
	  material_interactor mi(rp_,material,random_,scattering_optical_depth,
				 sampling_asymmetry_factor);
	  //std::cout << "mi set" << std::endl;
	  wall_interactor wi(nav_,rp_,random_,lag_);
	  //std::cout << "wi set" << std::endl;
	  //std::cout << mi.distance_to_scattering() << std::endl;
	  //std::cout << wi.distance_to_wall() << std::endl;
	  //std::cout << (mi.distance_to_scattering() < wi.distance_to_wall()) << std::endl;
	  if (mi.distance_to_scattering() < wi.distance_to_wall()) {
	    mi.deposite_energy_to_heat(mi.distance_to_scattering());
	    mi.scatter();
	    scattering_optical_depth = -log(random_(0,1)); //include likelihood?
	  } else {
	    //std::cout << "loop start" << std::endl;
	    mi.deposite_energy_to_heat(wi.distance_to_wall());
	    wi.interact_with_wall();
	    scattering_optical_depth -= material.scattering_optical_depth(wi.distance_to_wall());
	  }
	}
      }
    }
  };
}
}

#endif
