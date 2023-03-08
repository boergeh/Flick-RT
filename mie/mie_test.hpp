#include "mie.hpp"

namespace flick {
  begin_test_case(mie_test) {
    parameterized_monodisperesed_mie m({1,0},{1.33,0},500e-9);
    double r = 10e-6;
    m.set_radius(r);
    check_small(m.absorption_cross_section(), 1e-12);
    check_close(m.extinction_cross_section(), 2*3.14159*pow(r,2), 0.001);
    check(m.scattering_function(0,0)[0] > m.scattering_function(0,0)[1]);
  } end_test_case()
}
