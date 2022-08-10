#ifndef flick_mc_basic
#define flick_mc_basic

#include "../geometry/volume.hpp"
#include "../component/content.hpp"
#include "../component/emitter.hpp"
//#include "../coating/coating.hpp"

namespace flick {
  using half_infinite_box = geometry::half_infinite_box<content>;

  geometry::volume<content> world() {
    half_infinite_box space;
    half_infinite_box atmosphere;
    half_infinite_box bottom;

    space.name("space");
    atmosphere.name("atmosphere");
    bottom.name("bottom");

    atmosphere.content().receiver().activate();
    //atmosphere.content().receiver().outward_accepting();
    bottom.content().receiver().activate();
    bottom.content().set_coating<coating::grey_lambert>(1,0);
    
    space.move_by({0,0,10});
    atmosphere.move_by({0,0,1});
    atmosphere.insert(bottom);
    space.insert(atmosphere);
    
    return space;
  }
  class wall_interactor {
    geometry::navigator<content>* nav_;
    std::optional<pose> next_wall_intersection_;
    geometry::volume<content>* next_volume_;
    const coating::base* coating_;
    radiation_package* rp_;
    uniform_random rnd_;
    coating::lambert_angle_generator lag_;
    bool is_reflected_{false};
    bool is_transmitted_{true};
  public:
    wall_interactor(geometry::navigator<content>& nav,
		    radiation_package& rp) {
      nav_ = &nav;
      rp_ = &rp;
      next_wall_intersection_ = nav_->next_intersection(rp_->pose());
      next_volume_ = &nav_->next_volume(rp_->pose());
      auto& v = *next_volume_;
      if (is_moving_outward())
	v = nav_->current_volume();
      if (v.content().has_coating()) {
	coating_ = &v.content().coating();
	lag_ = coating::lambert_angle_generator{rnd_(0,1),rnd_(0,1)};
	is_reflected_=
	  (coating_!=nullptr && rnd_(0,1) < coating_->reflectivity());
	is_transmitted_=
	  (coating_==nullptr || rnd_(0,1) < coating_->transmissivity());
      }
    }
    bool will_hit() const {
      if (next_wall_intersection_.has_value())
	return true;
      return false;
    }
    void move_to_hit_point() const {
      rp_->move_to((*next_wall_intersection_).position());
    }    
    void interact() {
      unit_vector n = facing_surface_normal();
      if (is_reflected_) {
	change_orientation(n);
	change_radiation();
      } else if (is_transmitted_) {
	change_orientation(-n);
	nav_->go_to(*next_volume_);
	auto& re = nav_->current_volume().content().receiver();
	re.receive(*rp_, is_moving_inward());		
      } else {
	absorb_radiation_package();
      }      
    }
    auto& current_volume() const {
      return nav_->current_volume();
    }
    unit_vector facing_surface_normal() const {
      unit_vector n = (*next_wall_intersection_).z_direction();
      if (is_moving_outward())
	n = -n;
      return n;
    }
    void absorb_radiation_package() {
      rp_->scale_intensity(0);
    }
    bool is_moving_outward() const {
      if (nav_->is_moving_inward(*next_wall_intersection_,rp_->pose()))
	return false;
      return true;
    }    
    bool is_moving_inward() const {
      return !is_moving_outward();
    }    
    void change_orientation(const unit_vector& facing_surface_normal) {
      rp_->x_direction_parallel_with_plane_of_incidence(facing_surface_normal);
      rp_->traveling_direction(facing_surface_normal);
      rp_->rotate_about_local_z(lag_.azimuth_angle());
      rp_->rotate_about_local_x(lag_.polar_angle());
    }  
    void change_radiation() {
      double f = 1;
      if (is_reflected_) {
	double brdf = coating_->reflection_mueller_matrix().value(0,0);
	f = coating_->reflectivity()*brdf/lag_.brdf();
      } else if (is_transmitted_) {
	double brdf = coating_->transmission_mueller_matrix().value(0,0);
	f = coating_->transmissivity()*brdf/lag_.brdf();
      }
      rp_->scale_intensity(f);
    }
    radiation_package& leaving_radiation_package() {
      return *rp_;
    }
  };
  
  class mc_basic {
    geometry::volume<content>* outer_volume_;
    geometry::volume<content>* emitter_volume_;
    geometry::navigator<content> nav_;
    //emitter emitter_;
    radiation_package rp_;   
    uniform_random rnd_;
    
    bool lost_in_space() {
      if (nav_.next_volume(rp_.pose()).has_outer_volume())
      	return false;
      return true;
    }
  public:
    mc_basic(geometry::volume<content>& v) {
      outer_volume_ = &v;
      nav_ = geometry::navigator<content>(*outer_volume_);
      emitter_volume_ = &nav_.find("space");
    }
    void run() {
      emitter emitter{{0,0,5},stokes{1,0,0,0},3};
      emitter.set_wavelength<monocromatic>(500e-9);
      emitter.set_direction<unidirectional>(unit_vector{0,0,-1});
      while (!emitter.is_empty()) {
	nav_.go_to(*emitter_volume_);
	rp_ = emitter.emit();
	while (!rp_.is_empty() && !lost_in_space()) {
	  wall_interactor wi(nav_,rp_);
	  //particle_interactor pi(nav_,rp_);
	  if (wi.will_hit()) {	   
	    wi.move_to_hit_point();	   	  
	    wi.interact();         
	    rp_ = wi.leaving_radiation_package();
	    //nav_.go_to(wi.current_volume());
	  } //else if (pi.will_scatter()) {
	  // pi.move_to_scattering_event()
	    //pi.scatter()
	  //}
	  //else
	  //  nav_.go_to(nav_.next_volume(rp_.pose()));  
	}
      }
      auto& v = nav_.find("atmosphere");
      std::cout <<'\n'<< v.content().receiver(); 
    }
  };
}

#endif
