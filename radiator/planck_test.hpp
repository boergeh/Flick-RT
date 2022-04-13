#include "planck.hpp"

namespace flick {
  begin_test_case(planck_test) {
    double p = planck(5800).value(500e-9);
    check_close(p,8.5e13,10);
    
  } end_test_case()
}
