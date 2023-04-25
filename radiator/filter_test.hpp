#include "filter.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  begin_test_case(filter_test) {
    filter::gaussian g{440e-9,10e-9};
    pl_function f = read<pl_function>("./toa_solar.txt");
    f = transmit(f,filter::cut_ends(280e-9,400e-9));
    check_close(uv_index(f),uva_index(f)+uvb_index(f),0.3);
  } end_test_case()

}
