#include "profile.hpp"

namespace flick {
  begin_test_case(profile_test) {

    auto cp_o3 = concentration_profile("o3.txt");
    auto cp_co2 = concentration_profile("co2.txt");
    auto T = temperature_profile("temperature.txt");
    auto P = pressure_profile("pressure.txt");
    check_close(cp_o3.stp_thickness(),0.00348, 1);
    check_close(cp_co2.stp_thickness(),2.7, 1);
    check_close(P.value(5e3)/P.value(0),0.53, 1);
    
  } end_test_case()
}
