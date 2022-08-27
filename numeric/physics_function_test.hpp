#include "physics_function.hpp"
namespace flick {
  begin_test_case(physics_function_test) {
    henyey_greenstein hg{0};
    double P_iso = 1/(4*constants::pi);
    check_close(hg.phase_function(0),P_iso,1e-12,"a");
    check_small(hg.inverted_accumulated_angle(0),1e-12,"b");
    check_close(hg.inverted_accumulated_angle(1),constants::pi,1e-12,"c");
    hg = henyey_greenstein{0.9};
    check_small(hg.inverted_accumulated_angle(0),1e-7,"d");
    check_close(hg.inverted_accumulated_angle(1),constants::pi,1e-6,"e");
  } end_test_case()
}
