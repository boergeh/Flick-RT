#include "atmospheric_state.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  using namespace units;
  begin_test_case(atmospheric_state_test) {
    atmospheric_state state_1(296_K, 1013_hPa, 1);
    atmospheric_state state(296_K, 1013_hPa, 50);
    check_close(state_1.stp_thickness("o3"), state.stp_thickness("o3"));
    check_close(state_1.stp_thickness("o2"), state.stp_thickness("o2"));
    
    check_close(state.temperature(0_m), 296_K);
    check_close(state.stp_thickness("o3"), 3.48_mm, 1_pct);
    state.scale_to_stp_thickness("o3", 2.50_mm);
    check_close(state.stp_thickness("o3"), 2.50_mm);
    check_small(state.pressure(1000_km));
    
  } end_test_case()
}
