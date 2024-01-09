#include "pure_water.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(pure_water_test_A) {
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
    check_close(pw.real_refractive_index(),1.37,0.2_pct);
    pw.set_wavelength(1001_nm);
    check_close(pw.real_refractive_index(),1.326,0.2_pct);
    pw.set_wavelength(2000_nm);
    check_close(pw.real_refractive_index(),1.31,0.2_pct);
  } end_test_case()
  
   begin_test_case(pure_water_test_B) {
    // Zhang, X., and L. Hu. 2009. Estimating scattering of pure water
    // from density fluctuation of the refractive
    // index. Opt. Expr. 17: 1671- 1678.
    material::pure_water pw;
    unit_vector v90 = {0,1,0};   
    unit_vector v0 = {0,0,1};   
    pw.set_direction(v0);
    pw.set_wavelength(405e-9);
    double p90 = pw.mueller_matrix(v90).value(0,0);
    double scat_coef = pw.scattering_coefficient();
    //check_close(scat_coef, pow(129.0 / (pw.wavelength() * 1e9), 4.32));
    //check_close(p90,1/(4*3.14));
    //double vol_scat_90 = scat_coef*p90;
    //check_close(vol_scat_90, 2.90e-4, 2_pct);
  } end_test_case()
}
