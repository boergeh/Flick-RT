#include "air.hpp"

namespace flick {
  begin_test_case(air_test) {
    atmospheric_state state(300,1000e2);
    material::air a(state);
    //std::cout << a.scattering_optical_depth(100e3);
    a.set_wavelength(800e-9);
    //std::cout << a.scattering_optical_depth(100e3);
  } end_test_case()
}
