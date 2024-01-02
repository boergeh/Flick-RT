#include "time_point.hpp"
#include "../environment/unit_test.hpp"

namespace flick {
  begin_test_case(time_point_test) {
    time_point J0 = {-4713,11,24,12,0,0};
    check_small(J0.julian_date());
    check_small(time_point().julian_date());
    time_point  J2000 = {2000,1,1,12,0,0};
    check_small(J2000.julian_date()-2451545.0);
    //std::cout << J2000;
  } end_test_case()
}
