#include "planck.hpp"
#include "toa_solar.hpp"

namespace flick {
  begin_test_case(planck_test) {
    double T = 5800;
    auto p = planck(T);
    check_close(p.irradiance(500e-9),8.5e13,10.0_pct);
    check_close(p.inverted_accumulated_wavelength(0.1),400e-9,6.0_pct);
    auto irr = p.irradiance_spectrum(150);
    check_close(irr.integral(),constants::sigma*pow(T,4),1.0_pct);
  } end_test_case()

  begin_test_case(toa_solar_test) {
    radiator::toa_solar s;
    check_close(s.spectrum(170).integral(),1361,1.0_pct);
  } end_test_case()
}
