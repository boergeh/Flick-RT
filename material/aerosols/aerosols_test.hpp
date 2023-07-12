#include "../../environment/unit_test.hpp"
#include "aerosols.hpp"

namespace flick {
  using namespace flick;
  begin_test_case(aerosols_test) {
    material::rural_aerosols ra;
    ra.set(pose{{0,0,1},unit_vector{0,0,1}});
    ra.set(wavelength{500e-9});
    check_close(ra.absorption_coefficient(),1e-5,40.0_pct);
  } end_test_case()
}
