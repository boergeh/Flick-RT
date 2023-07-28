#include "profile.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(profile_test) {
    using namespace units;
    auto cp_o3 = concentration_profile("o3.txt");
    auto cp_co2 = concentration_profile("co2.txt");
    auto T = temperature_profile("temperature.txt");
    auto P = pressure_profile("pressure.txt");
    check_close(T.value(0_m), 288.2_K);
    check_close(cp_o3.stp_thickness(), 3.48_mm, 1_pct);
    check_close(cp_co2.stp_thickness(), 2.7_m, 1_pct);
    check_close(P.value(5_km)/P.value(0_m), 0.53, 1_pct);
  } end_test_case()
}
