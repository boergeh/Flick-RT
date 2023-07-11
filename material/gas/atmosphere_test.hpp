#include "atmosphere.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    o3_cross_section c;
    atmosphere_state state(300,1000e2);
    state.remove_gas("h2o");
  } end_test_case()
}
