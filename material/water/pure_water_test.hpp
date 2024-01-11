#include "pure_water.hpp"
#include "../../numeric/units.hpp"
#include "../../numeric/constants.hpp"

namespace flick {
  begin_test_case(pure_water_test_A) {
    // Zhang, X., and L. Hu. 2009. Estimating scattering of pure water
    // from density fluctuation of the refractive
    // index. Opt. Expr. 17: 1671- 1678, Table 1.
    double S = 0;
    double wl = 366e-9;
    double PMH = 4.500e-4;
    double T = constants::to_kelvin(20);
    material::zhang_hu_scattering zhs(S,T,wl);
    check_close(zhs.volume_scattering_90(), PMH, 0.1_pct);
    wl = 578e-9;
    PMH = 0.650e-4;
    material::zhang_hu_scattering zhs2(S,T,wl);
    check_close(zhs2.volume_scattering_90(), PMH, 0.2_pct);
    
    S = 0;
    wl = 300e-9;
    T = constants::to_kelvin(20);
    material::zhang_hu_scattering zhs3(S,T,wl);
    //double b_morel = pow(129.0 / (wl * 1e9), 4.32);
    // https://www.oceanopticsbook.info/view/optical-constituents-of-the-ocean/water
    check_close(zhs3.scattering_coefficient(), 0.018, 0.5_pct);
    S = 35;
    T = constants::to_kelvin(0);
    material::zhang_hu_scattering zhs4(S,T,wl);
    check_close(zhs4.scattering_coefficient(), 0.025, 0.5_pct);
    
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
    check_close(pw.real_refractive_index(),1.37,0.2_pct);
    pw.set_wavelength(1001_nm);
    check_close(pw.real_refractive_index(),1.326,0.2_pct);
    pw.set_wavelength(2000_nm);
    check_close(pw.real_refractive_index(),1.31,0.2_pct);
  } end_test_case()
  
   begin_test_case(pure_water_test_C) {
    //return pow(129.0 / (wavelength() * 1e9), 4.32) * volume_fraction_;
    material::pure_water pw;
    unit_vector v0 = {0,0,1};   
    unit_vector v90 = {0,1,0};   
    pw.set_direction(v0);
    double p0 = pw.mueller_matrix(v0).value(0,0);
    double p90 = pw.mueller_matrix(v90).value(0,0);
    check_close(p90/p0,0.5,4_pct);
  } end_test_case()
}
