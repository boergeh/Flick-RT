#include "physics_function.hpp"

namespace flick {
  begin_test_case(physics_function_test) {
    using namespace constants;
    henyey_greenstein hg{0};
    double P_iso = 1/(4*pi);
    check_close(hg.phase_function(0), P_iso);
    check_small(hg.inverted_accumulated_angle(0));
    check_close(hg.inverted_accumulated_angle(1), pi);
    hg = henyey_greenstein{0.9};
    check_small(hg.inverted_accumulated_angle(0),1e-7);
    check_close(hg.inverted_accumulated_angle(1), pi,1e-6_pct);
  } end_test_case()
}
