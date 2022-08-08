#ifndef flick_mc_basic
#define flick_mc_basic

#include "../geometry/volume.hpp"
#include "../component/content.hpp"
#include "../component/emitter.hpp"
//#include "../coating/coating.hpp"

namespace flick {
  using half_infinite_box = geometry::half_infinite_box<content>;
  void mc_basic() {
    half_infinite_box space;
    half_infinite_box atmosphere;
    half_infinite_box bottom;

    space.name("space");
    atmosphere.name("atmosphere");
    bottom.name("bottom");

    bottom.content().receiver().activate();
    //bottom.content().set_coating<coating::white_lambert>();
    //bottom.content().set_coating<coating::white_lambert>();
    //space.content().set_coating<coating::grey_lambert>(1,0);
    //atmosphere.content().set_coating<coating::grey_lambert>(1,0);
    bottom.content().set_coating<coating::grey_lambert>(1,0);

    //coating::white_lambert wla;
    //std::cout << wla.reflectivity();
    
    space.move_by({0,0,10});
    atmosphere.move_by({0,0,1});
    atmosphere.insert(bottom);
    space.insert(atmosphere);

    emitter e{{0,0,1},stokes{1,0,0,0},3};
    e.set_wavelength<monocromatic>(500e-9);
    e.set_direction<unidirectional>(unit_vector{0,0,-1});

    auto nav = geometry::navigator<content>(space);
    auto& atm = nav.go_to("atmosphere");
    
    std::cout << std::endl;
    uniform_random rnd_;
    while (!e.is_empty()) {
      nav.go_to(atm);
      radiation_package rp = e.emit();
      while (!rp.is_empty() && nav.current_volume().has_outer_volume()) {
	auto& cv = nav.current_volume();
	auto& nv = nav.next_volume(rp.pose());
	auto ni = nav.next_intersection(rp.pose());

	if (ni.has_value())
	  rp.move_to((*ni).position());

	std::cout << "current volume: "<<cv.name() << std::endl;
	std::cout << "position: "<<rp.pose().position().z() << std::endl;
	std::cout << "direction: "<<rp.pose().direction() << std::endl;
		

	if (ni.has_value() && nav.is_entering(*ni,rp.pose())
	    && nv.content().has_coating()) {
	  // for reflection
	  auto& c = nv.content().coating();
	  
	  unit_vector surface_normal = (*ni).z_direction();
	  std::cout << "surfnorm " << surface_normal << std::endl;
	  coating::lambert_angle_generator lag;
	  //coating::lambert_angle_generator lag(rnd(0,1),rnd(0,1));
	  //rp.x_direction_parallel_with_plane_of_incidence(surface_normal);
	  /*
	  rp.traveling_direction(surface_normal);
	  rp.rotate_about_local_z(lag.azimuth_angle());
	  rp.rotate_about_local_x(lag.polar_angle());
	  double brdf = c.reflection_mueller_matrix().value(0,0);
	  double f = c.reflectivity()*brdf/lag.brdf();
	  rp.scale_intensity(f);
	  nv = nav.next_volume(rp.pose());
	  */	  
//rp.traveling_direction(dg_.lambertian(relative_normal()));
	  //rp.change_polarization(mueller,interaction_plane_angle)
	
	}
	/*
	if (cc. > rnd_(0,1)) { 
	} else { 
	  rp.scale_intensity(transmittance_);
	  rp.traveling_direction(dg_.lambertian(-relative_normal()));
	}
	*/
	nav.go_to(nv);
	nv.content().receiver().receive(rp);
	std::cout << std::endl;
      }
    }
    auto bot = nav.go_to("bottom");
    std::cout<<bot.content().receiver().received_packages();
    std::cout<<bot.content().receiver();
    
  }
}

#endif
