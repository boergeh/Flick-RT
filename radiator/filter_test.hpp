#include "filter.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  begin_test_case(filter_test) {
    double p = 0.1;
    check_close(filter::erythema().transmittance(298e-9),1,p,"a");
    check_close(filter::erythema().transmittance(310e-9),0.07447,p,"b");
    check_close(filter::erythema().transmittance(350e-9),0.000707946,p,"c");
    pl_function f = read<pl_function>("./toa_solar.txt");
    auto f2 = transmit(f,filter::cut_ends(280e-9,400e-9));
    check_close(uv_index(f2),uva_index(f2)+uvb_index(f2),0.3,"d");
    double wl0 = 0.5e-6;
    double fwhm = 5e-9;
    check_close(gaussian_mean(f,wl0,fwhm),square_mean(f,wl0,fwhm),1,"e");
  } end_test_case()

}
