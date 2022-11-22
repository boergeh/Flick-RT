#include "atmosphere_state.hpp"

namespace flick {
  begin_test_case(atmosphere_state_test) {
    atmosphere_state a(300,1000*1e2);
    check_close(a.temperature(0),300,1e-12);
    check_close(a.stp_thickness(gas::o3),344e-5,1);
    check_close(a.stp_thickness(gas::air),8000,1);
    check_close(a.temperature(100e3),203,1);
    check_small(a.pressure(200e3),1e-3);
    
  } end_test_case()
}
