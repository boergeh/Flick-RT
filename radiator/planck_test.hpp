#include "planck.hpp"

namespace flick {
  begin_test_case(planck_test) {
    double T = 5800;
    auto p = planck(T);
    check_close(p.irradiance(500e-9),8.5e13,10);
    check_close(p.density_wavelength(0.1),400e-9,6);
    auto irr = p.irradiance_spectrum(25);
    check_close(irr.integral(),constants::sigma*pow(T,4),5);
  } end_test_case()
}
