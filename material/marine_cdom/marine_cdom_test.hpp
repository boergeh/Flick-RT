#include "marine_cdom.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(marine_cdom_test) {
    using namespace units;
    material::marine_cdom mcdom("HF22_D001",1);
    mcdom.set_wavelength(550e-9);
    check_close(mcdom.absorption_coefficient(),0.001,0.1_pct);
    check_small(mcdom.scattering_coefficient());
  } end_test_case()
}
