#include "filter.hpp"

namespace flick {
  begin_test_case(filter_test) {
    filter::gaussian g{440e-9,10e-9,100};
    check_close(g.spectrum(100).integral(),1,1);
  } end_test_case()

}
