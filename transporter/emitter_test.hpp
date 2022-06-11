#include "emitter.hpp"
namespace flick {
  begin_test_case(emitter_test) {
    emitter em({0,0,0},stokes{1,0,0,0},1e3);
    //em.placement({{1,1,1},{0,0}});
    em.set_wavelength<monocromatic>(500e-9);
    em.set_direction<unidirectional>(unit_vector{0,0});
    em.set_direction<conic>(constants::pi/2, unit_vector{0,0,-1});
    check(em.emit().pose().direction().z() <= 0);    
    //std::cout<< em;
  } end_test_case()
}
