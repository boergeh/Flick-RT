#include "pure_water.hpp"
#include "../../numeric/units.hpp"
#include "../../numeric/constants.hpp"

namespace flick {
  begin_test_case(pure_water_test_A) {
    // Zhang, X., and L. Hu. 2009. Estimating scattering of pure water
    // from density fluctuation of the refractive
    // index. Opt. Expr. 17: 1671- 1678, Table 1.
    double S = 0;
    double T = constants::to_kelvin(20);
    material::water::scattering s1(S,T);
    check_close(s1.vsf90_density(366e-9), 4.500e-4, 0.3_pct);
    material::water::scattering s2(S,T);
    double bench = 0.650e-4;
    double wl = 578e-9;
    check_close(s2.vsf90_density(wl), bench, 0.2_pct);
      
    // https://www.oceanopticsbook.info/view/
    // optical-constituents-of-the-ocean/water
    S = 0;
    T = constants::to_kelvin(20);
    bench = 0.018;
    wl = 300e-9;
    material::water::scattering s3(S,T);
    check_close(s3.coefficient(wl), bench, 0.3_pct);

    S = 35;
    T = constants::to_kelvin(20);
    wl = 300e-9;
    material::water::scattering s4(S,T);
    double flick = s4.coefficient(wl);
    double oo_web_book = 0.0235;
    double zhang_hu_matlab = 0.02361536;
    double morel = pow(129.0 / (wl * 1e9), 4.32);
    check_close(flick, oo_web_book, 1_pct,
		"Figure in ocean optics web book");
    check_close(flick, zhang_hu_matlab, 0.0001_pct,
		"Zhang and Hu (2009) Matlab script");
    check_close(flick, morel, 10_pct,
		"Morel parameterization");
  } end_test_case()
  
  begin_test_case(pure_water_test_B) {  
    using namespace units;
    material::pure_water pw;
    check_fast(70*cpu_duration());
    pw.set_wavelength(500_nm);
    check_close(pw.absorption_coefficient(),0.02,5.0_pct);
    pw.salinity(pl_function{35_psu});
    pw.temperature({273_K});
    check_close(pw.absorption_coefficient(),0.021,2.0_pct);
    pw.set_wavelength(835_nm);
    check_close(pw.absorption_coefficient(),2.99,1.0_pct);
    pw.salinity(pl_function{0});
    check_close(pw.absorption_coefficient(),3.02,1.0_pct);
    pw.temperature({273+30_K});
    check_close(pw.absorption_coefficient(),3.47,1.0_pct);

    pw.set_wavelength(270_nm);
    pw.temperature({273_K});
    pw.salinity(0_psu);
    check_close(pw.real_refractive_index(),1.37,0.5_pct);
    pw.set_wavelength(1001_nm);
    check_close(pw.real_refractive_index(),1.326,0.2_pct);
    pw.set_wavelength(2000_nm);
    check_close(pw.real_refractive_index(),1.31,0.2_pct);
  } end_test_case()
  
  begin_test_case(pure_water_test_C) {
    using namespace units;
    material::pure_water pw;
    pw.temperature({273.15_K});
    pw.salinity(0_psu);
    pw.set_wavelength(340e-9);
    check(pw.absorption_coefficient() > 0);
  } end_test_case()
  
  begin_test_case(pure_water_test_D) {
    // Assert small noise in uv
    using namespace units;
    material::pure_water pw;
    pw.temperature({273.15_K});
    pw.salinity(0_psu);
    size_t n = 100;
    auto wl = range(340e-9,350e-9,n).linspace();
    //auto wl = range(540e-9,550e-9,n).linspace();
    double y = 0;
    for (size_t i=0; i<wl.size(); i++) {
      pw.set_wavelength(wl[i]);
      y += pw.absorption_coefficient();
      //std::cout << wl[i] << " " << y << std::endl; 
    }
    double mean = y/n;
    y = 0;
    for (size_t i=0; i<wl.size(); i++) {
      pw.set_wavelength(wl[i]);
      y += fabs(pw.absorption_coefficient()-mean);
    }
    double variation = y/n;
    check_small(variation/mean, 0.1);
  } end_test_case()
  
   begin_test_case(pure_water_test_E) {
    double S = 35;
    double T = constants::to_kelvin(0);
    material::pure_water pw(S,T);
    unit_vector v0 = {0,0,1};   
    unit_vector v90 = {0,1,0};   
    pw.set_direction(v0);
    double p0 = pw.mueller_matrix(v0).value(0,0);
    double p90 = pw.mueller_matrix(v90).value(0,0);
    check_close(p90/p0,0.5,4_pct);
    pw.set_wavelength(400e-9);
    check_close(p0*pw.scattering_coefficient(),8.2053e-4,0.001_pct,"Zhang-Hu Matlab script");
    check_close(p90*pw.scattering_coefficient(),4.2626e-4,0.001_pct,"Zhang-Hu Matlab script");
  } end_test_case()
}
