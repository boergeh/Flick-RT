#include "range.hpp"

namespace flick {
  begin_test_case(range_test) {
    check_small(range(2,3,10).logspace()[8]
		- valarray_range(2,3,10).logspace()[8]);
  } end_test_case()
}
