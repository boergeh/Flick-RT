#include "coating.hpp"

namespace flick {  
  using namespace constants;
  begin_test_case(coating_test) {
    coating::white_lambert l;
    //check_close(l.reflection_polar_angle(), pi/2, 1e-12);
    //coating::fresnel f({1.33,0});
    //check_small(f.transmission_polar_angle(),1e-12);
    //check_small(f.reflection_polar_angle(),1e-12);
    //auto m = f.reflection_mueller_matrix();
    //std::cout << m;
    //check(m(0).value > 0);
  } end_test_case()
}
