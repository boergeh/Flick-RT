#ifndef flick_polarization_stokes_test
#define flick_polarization_stokes_test

#include "stokes.hpp"

namespace flick {
  begin_test_case(stokes_test) {
    stokes s{{1,-1,0,0}};
    check_close(s.Q(),-1,1e-12);
  } end_test_case()
}

#endif
