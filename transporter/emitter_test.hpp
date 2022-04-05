#include "emitter.hpp"
namespace flick {
  begin_test_case(emitter_test) {
    emitter em(1000,stokes{1,0,0,0});
    //em.placement({{1,1,1},{0,0}});
    em.wavelength<monocromatic>(500e-9);
    em.direction<unidirectional>(unit_vector{0,0});
    em.direction<conic>(pi/2, unit_vector{0,0,-1});
    check(em.emit().pose().direction().z() <= 0);    
    //std::cout<< em;
  } end_test_case()
}
