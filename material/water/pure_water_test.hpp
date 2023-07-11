#include "pure_water.hpp"

namespace flick {
  using namespace flick;
  begin_test_case(pure_water_test) {
    material::pure_water pw;
    pw.set(wavelength{500e-9});
    check_close(pw.absorption_coefficient(),0.02,5);
    pw.salinity(pl_function{35});
    pw.temperature({273});
    check_close(pw.absorption_coefficient(),0.021,2);
    pw.set(wavelength{835e-9});
    check_close(pw.absorption_coefficient(),2.99,1);
    pw.salinity(pl_function{0});
    check_close(pw.absorption_coefficient(),3.02,1);
    pw.temperature({273+30});
    check_close(pw.absorption_coefficient(),3.47,1);

    pw.set(wavelength{270e-9});
    pw.temperature({273+0});
    pw.salinity(0);
    check_close(pw.real_refractive_index(),1.37,0.2);
    pw.set(wavelength{1001e-9});
    check_close(pw.real_refractive_index(),1.326,0.2);
    pw.set(wavelength{2000e-9});
    check_close(pw.real_refractive_index(),1.31,0.2);
  } end_test_case()
}
