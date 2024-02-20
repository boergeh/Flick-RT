#include "marine_cdom.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(marine_cdom_test) {
    using namespace units;
    material::marine_cdom mcdom("ECOSENS_HF22_D1",1);
    mcdom.set_wavelength(301e-9);
    check_close(mcdom.absorption_coefficient(),1.724295,0.001_pct);
    check_small(mcdom.scattering_coefficient());
  } end_test_case()
}
