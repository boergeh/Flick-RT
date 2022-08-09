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
    bottom.content().receiver().activate();
    bottom.content().set_coating<coating::grey_lambert>(1,0);
    
    space.move_by({0,0,10});
    atmosphere.move_by({0,0,1});
    atmosphere.insert(bottom);
    space.insert(atmosphere);
    
    return space;
  }
  class mc_basicm {
    geometry::volume<content>* outer_volume_;
    geometry::volume<content>* emitter_volume_;
    geometry::navigator<content> nav_;
    emitter emitter_;
    radiation_package rp_;
    std::optional<pose> next_wall_intersection_;
    const coating::base* coating_;
    uniform_random rnd_;
    
    bool lost_in_space() {
      if (nav_.next_volume(rp_.pose()).has_outer_volume())
      	return false;
      return true;
    }
    //    void check_for_wall_intersection() {
    //  next_wall_intersection_ = nav_.next_intersection(rp_.pose());
    // }
    bool will_hit_wall() {
      next_wall_intersection_ = nav_.next_intersection(rp_.pose());
      if (next_wall_intersection_.has_value())
	return true;
      return false;
    }
    void move_to_wall() {
      rp_.move_to((*next_wall_intersection_).position());
    }
    void interact_with_wall() {
      auto& v = nav_.next_volume(rp_.pose());
      if (moving_outward())
	v = nav_.current_volume(); 	
      if (v.content().has_coating())
	coating_ = &v.content().coating();
      if (is_reflected()) {
	change_orientation();
	//change_polarization();
	//} else if (is_transmitted())
      }
      
    }
    void change_orientation() {
      unit_vector surface_normal = (*next_wall_intersection_).z_direction();
      if (moving_outward())
	surface_normal = -surface_normal;
      coating::lambert_angle_generator lag(rnd_(0,1),rnd_(0,1));
      rp_.x_direction_parallel_with_plane_of_incidence(surface_normal);	  
      rp_.traveling_direction(surface_normal);
      rp_.rotate_about_local_z(lag.azimuth_angle());
      rp_.rotate_about_local_x(lag.polar_angle());
    }

    bool moving_outward() {
      if (nav_.is_moving_inward(*next_wall_intersection_,rp_.pose()))
	return false;
      return true;
    }
    
    bool is_reflected() {
      if (rnd_(0,1) < coating_->reflectivity())
	return true;
      return false;
    }
  public:
    mc_basicm(geometry::volume<content>& v) {
      outer_volume_ = &v;
      nav_ = geometry::navigator<content>(*outer_volume_);
      emitter_volume_ = &nav_.find("space");
    }
    void run() {
      while (!emitter_.is_empty()) {
	nav_.go_to(*emitter_volume_);
	rp_ = emitter_.emit();
	while (!rp_.is_empty() && !lost_in_space()) {
	  if (will_hit_wall()) {
	    move_to_wall();	   
	    interact_with_wall();
	  //}  else if (will_scatter) {
	  // move_to_scattering_event()
	    //scatter()
	  }
	}
      }
      /*
      while (!emitter.is_empty()) {
	nav.go_to(space_r);
	radiation_package rp = emitter.emit();
	while (!rp.is_empty() && nav.next_volume(rp.pose()).has_outer_volume()) {
	  auto& cv = nav.current_volume();
	  auto& nv = nav.next_volume(rp.pose());
	  auto ni = nav.next_intersection(rp.pose());
	  if (ni.has_value())
	    rp.move_to((*ni).position());
	  if (ni.has_value() && nav.is_entering(*ni,rp.pose())
	      && nv.content().has_coating()) {
	    auto& c = nv.content().coating();
	  
	    unit_vector surface_normal = (*ni).z_direction();
	    coating::lambert_angle_generator lag;
	    //coating::lambert_angle_generator lag(rnd(0,1),rnd(0,1));

	    rp.x_direction_parallel_with_plane_of_incidence(surface_normal);	  
	    rp.traveling_direction(surface_normal);
	    rp.rotate_about_local_z(lag.azimuth_angle());
	    rp.rotate_about_local_x(lag.polar_angle());
	    
	    double brdf = c.reflection_mueller_matrix().value(0,0);
	    double f = c.reflectivity()*brdf/lag.brdf();
	    rp.scale_intensity(f);
	    
	    nv = nav.next_volume(rp.pose());
	  }
	
	  nav.go_to(nv);
	  nv.content().receiver().receive(rp);
	}
      }
      */
    }
  };

  void mc_basic() {
    geometry::volume<content> w = world();
    mc_basicm mc{w};
    
    /*
    emitter emitter{{0,0,5},stokes{1,0,0,0},3};
    emitter.set_wavelength<monocromatic>(500e-9);
    emitter.set_direction<unidirectional>(unit_vector{0,0,-1});

    auto w = world();
    auto nav = geometry::navigator<content>(w);
    auto& space_r = nav.find("space");
    
    uniform_random rnd;
    while (!emitter.is_empty()) {
      nav.go_to(space_r);
      radiation_package rp = emitter.emit();
      while (!rp.is_empty() && nav.next_volume(rp.pose()).has_outer_volume()) {
	auto& cv = nav.current_volume();
	auto& nv = nav.next_volume(rp.pose());
	auto ni = nav.next_intersection(rp.pose());
	if (ni.has_value())
	  rp.move_to((*ni).position());
	if (ni.has_value() && nav.is_entering(*ni,rp.pose())
	    && nv.content().has_coating()) {
	  // for reflection
	  auto& c = nv.content().coating();
	  
	  unit_vector surface_normal = (*ni).z_direction();
	  coating::lambert_angle_generator lag;
	  //coating::lambert_angle_generator lag(rnd(0,1),rnd(0,1));

	  rp.x_direction_parallel_with_plane_of_incidence(surface_normal);	  
	  rp.traveling_direction(surface_normal);
	  rp.rotate_about_local_z(lag.azimuth_angle());
	  rp.rotate_about_local_x(lag.polar_angle());

	  double brdf = c.reflection_mueller_matrix().value(0,0);
	  double f = c.reflectivity()*brdf/lag.brdf();
	  rp.scale_intensity(f);

	  nv = nav.next_volume(rp.pose());
	}
	
	nav.go_to(nv);
	nv.content().receiver().receive(rp);
      }
    }
    auto bot = nav.go_to("atmosphere");
    std::cout<<bot.content().receiver().received_packages();
    std::cout<<bot.content().receiver();
    */
  }
}

#endif
