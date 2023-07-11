#include "stokes.hpp"

namespace flick {
  begin_test_case(stokes_test) {
    stokes s{{1,-1,0,0}};
    check_close(s.Q(),-1);
  } end_test_case()
}

