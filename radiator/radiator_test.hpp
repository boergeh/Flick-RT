#include "planck.hpp"
#include "toa_solar.hpp"

namespace flick {
  begin_test_case(planck_test) {
    double T = 5800;
    auto p = radiator::planck(T);
    check_close(p.irradiance(500e-9),8.5e13,10);
    check_close(p.spectrum(100).integral(),constants::sigma*pow(T,4),1);
  } end_test_case()

  begin_test_case(toa_solar_test) {
    radiator::toa_solar s;
    check_close(s.spectrum(170).integral(),1361,1);
  } end_test_case()
}
