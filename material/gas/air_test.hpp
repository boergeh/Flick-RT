#include "air.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(air_test) {
    using namespace units;
    atmospheric_state state(290_K, 1019_hPa, 50);
    state.remove_gas("co2");
    state.remove_gas("h2o");
    state.remove_gas("o2");
    state.scale_to_stp_thickness("o3", 2.95_mm);

    material::air air(state);
    air.set_wavelength(310_nm);
    air.set_position({0,0,17_m});
    air.set_direction({0,0,1});

    // Van Weele M, et al., Journal of Geophysical Research:
    // Atmospheres, 105(D4), pp.4915-4925, Table 5.
    check_close(air.scattering_optical_depth(120_km), 1.059, 0.2_pct);
    check_close(air.absorption_optical_depth(120_km), 0.687, 0.2_pct);

  } end_test_case()
}
