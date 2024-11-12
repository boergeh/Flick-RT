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
  
  begin_test_case(filter_test_C) {
    check_close(filter::cone_L().transmittance(600e-9),0.83399,0.1_pct);
    check_close(filter::tristimulus_y().transmittance(550e-9),0.989023,0.01_pct);
    size_t n = 500;
    pp_function E(range(380e-9,800e-9,n).linspace(),std::vector<double>(n,1e14));
    for (size_t i = 0; i<3; i++) {
      check_close(tristimulus(E)[i], 1./3, 0.001_pct);
    }
    // https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D
    pp_function A = radiator::planck(2856).spectrum(2000);
    check_close(tristimulus(A)[0], 0.44758, 1.2_pct);
    check_close(tristimulus(A)[1], 0.40745, 1.2_pct);
    
    // https://www.oceanopticsbook.info/view/photometry-and-visibility/from-xyz-to-rgb
    check_close(rgb(E)[0],255./255,0.2_pct);
    check_close(rgb(E)[1],201./255,0.2_pct);
    check_close(rgb(E)[2],192./255,0.3_pct);
  } end_test_case()
}
