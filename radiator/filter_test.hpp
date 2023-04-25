#include "filter.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  begin_test_case(filter_test) {
    //filter::gaussian g{440e-9,10e-9};
    double p = 0.1;
    check_close(filter::erythema().transmittance(298e-9),1,p);
    check_close(filter::erythema().transmittance(310e-9),0.07447,p);
    check_close(filter::erythema().transmittance(350e-9),0.0006839,p);
    pl_function f = read<pl_function>("./toa_solar.txt");
    auto f2 = transmit(f,filter::cut_ends(280e-9,400e-9));
    check_close(uv_index(f2),uva_index(f2)+uvb_index(f2),0.3);
    check_close(gaussian_mean(f,1e-6,0.1e-6),f.value(1e-6),0.1);
  } end_test_case()

}
