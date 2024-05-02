#include "marine_cdom.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(marine_cdom_test) {
    using namespace units;
    material::marine_cdom mcdom("ECOSENS_HF22_D1",1);
    mcdom.set_wavelength(301e-9);
    check_close(mcdom.absorption_coefficient(),1.724);
    check_small(mcdom.scattering_coefficient());
    mcdom.set_wavelength(780e-9);
    check_close(mcdom.absorption_coefficient(),-0.0118, 0.1_pct);
    check_small(mcdom.scattering_coefficient());
  } end_test_case()
}
