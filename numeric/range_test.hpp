#include "range.hpp"
namespace flick {
  begin_test_case(range_test) {
    //std::cout << range(2,3,10).linspace();
    //std::cout << valarray_range(2,3,10).logspace();
    check_small(range(2,3,10).logspace()[8]-valarray_range(2,3,10).logspace()[8],1e-15);
  } end_test_case()
}
