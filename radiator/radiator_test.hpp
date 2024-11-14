#include "planck.hpp"
#include "toa_solar.hpp"
#include "cie_d65.hpp"

namespace flick {
  begin_test_case(planck_test) {
    double T = 5800;
    auto p = radiator::planck(T);
    check_close(p.irradiance(500e-9),8.5e13,10.0_pct);
    check_close(p.importance_spectrum(400).integral(),constants::sigma*pow(T,4),0.1_pct);
  } end_test_case()

  begin_test_case(toa_solar_test) {
    radiator::toa_solar s;
    check_close(s.importance_spectrum(170).integral(),1361,1.0_pct);
  } end_test_case()
  
  begin_test_case(cie_d65_test) {
    radiator::cie_d65 s;
    check_close(s.spectrum().integral(),40852.20400e-9,0.1_pct);
  } end_test_case()
}
