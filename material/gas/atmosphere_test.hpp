#include "atmosphere.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    atmosphere_state state(300,1000e2);
    material::atmosphere a(state);
    //std::cout << a.scattering_optical_depth(100e3);
    a.set_wavelength(800e-9);
    //std::cout << a.scattering_optical_depth(100e3);
  } end_test_case()
}
