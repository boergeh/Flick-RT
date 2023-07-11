#include "content.hpp"

namespace flick {
  begin_test_case(content_test) {
    emitter em({0,0,0},stokes{1,0,0,0},1e3);
  } end_test_case()
}
