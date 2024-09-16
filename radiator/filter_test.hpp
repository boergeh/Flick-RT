#include "filter.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  begin_test_case(filter_test_A) {
    double p = 0.1;
    check_close(filter::erythema().transmittance(298e-9),1,p);
    check_close(filter::erythema().transmittance(310e-9),0.07447,p);
    check_close(filter::erythema().transmittance(350e-9),0.000707946,p);
    pl_function f = read<pl_function>("./toa_solar.txt");
    auto f2 = transmit(f,filter::cut_ends(280e-9,400e-9));
    check_close(uv_index(f2),uva_index(f2)+uvb_index(f2),0.3_pct);
    double wl0 = 0.5e-6;
    double fwhm = 5e-9;
    check_close(gaussian_mean(f,wl0,fwhm),square_mean(f,wl0,fwhm),1.0_pct);
  } end_test_case()
  
  begin_test_case(filter_test_B) {
    check(filter::sentinel3(300e-9).closest_srf()==0);
    check(filter::sentinel3(400e-9).transmittance(400e-9) > 1e-5);
    check_small(filter::sentinel3(300e-9).transmittance(420e-9));
    check_small(filter::sentinel3(300e-9).transmittance(380e-9));
    check(filter::sentinel3(1090e-9).closest_srf()==20);
    check(filter::sentinel3(938e-9).closest_srf()==19);
  } end_test_case()
}
