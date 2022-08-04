#ifndef flick_mc_basic
#define flick_mc_basic

#include "../geometry/volume.hpp"
#include "../component/content.hpp"
#include "../component/emitter.hpp"

namespace flick {
  using half_infinite_box = geometry::half_infinite_box<content>;
  void mc_basic() {
    half_infinite_box space;
    half_infinite_box atmosphere;
    half_infinite_box bottom;

    bottom.content().activate_receiver();
    space.move_by({0,0,10});
    atmosphere.move_by({0,0,1});
    atmosphere.insert(bottom);
    space.insert(atmosphere);

    emitter e{{0,0,1},stokes{1,0,0,0},10};
    e.set_wavelength<monocromatic>(500e-9);
    e.set_direction<unidirectional>(unit_vector{0,0,-1});

    auto nav = geometry::navigator<content>(space);
    nav.go_inward();

    radiation_package rp = e.emit();

    while (!e.is_empty() && !rp.is_empty() &&
	   nav.current_volume().has_outer_volume()) {
	auto v = nav.next_volume(rp.pose());
	auto p = nav.next_intersection(rp.pose());
	//rp.move_to(p);
	//nav.go_to_next_volume(rp.pose());
	//std::cout << "pose "<<rp.pose() << std::endl;
	//std::cout << "rp "<< rp << std::endl;
	
	rp = e.emit();
    }
  };
}

#endif
