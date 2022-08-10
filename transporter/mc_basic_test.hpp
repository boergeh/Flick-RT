#include "mc_basic.hpp"

namespace flick {  
  begin_test_case(mc_basic_test) {
    geometry::volume<content> w = world();
    mc_basic mcb{w};
    mcb.run();
  } end_test_case()
}
