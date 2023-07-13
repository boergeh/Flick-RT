#include "pure_ice.hpp"

namespace flick {
  begin_test_case(pure_ice_test) {
    using namespace flick;
    material::pure_ice ice;
    ice.set_wavelength(400e-9);
    check_close(ice.absorption_coefficient(),7.5e-4,5.0_pct);
    check_close(ice.real_refractive_index(),1.31,5.0_pct);
  } end_test_case()
}
