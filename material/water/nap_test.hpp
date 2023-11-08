#include "nap.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(nap_test) {
    using namespace units;
    material::nap p;
    p.set_wavelength(550e-9);
    p.mass_concentration(1e-3);
    check_close(p.scattering_coefficient(),0.5,30);
    p.set_wavelength(443e-9);
    check_close(p.absorption_coefficient(),0.03,30);
    auto p_forward = p.mueller_matrix(unit_vector{0.001,0}).value(0,0);
    auto unit_back = unit_vector{constants::pi,0};
    auto p_backward = p.mueller_matrix(unit_back).value(0,0);
    check(p_forward/p_backward > 1e4);
  } end_test_case()
}
