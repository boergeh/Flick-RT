#include "atmosphere_state.hpp"

namespace flick {
  begin_test_case(atmosphere_state_test) {
    atmosphere_state a(300,1000*1e2,3);
    check_close(a.temperature(0),300);
    check_close(a.stp_thickness("o3"),348e-5,1);
    check_small(a.pressure(200e3),1e-3);
    a.remove_gas("h2o");
  } end_test_case()
}
