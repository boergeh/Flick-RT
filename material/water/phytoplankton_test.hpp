#include "phytoplankton.hpp"
#include "../../numeric/units.hpp"
//#include "../material.hpp"

namespace flick {
  begin_test_case(phytoplankton_test) {
    using namespace units;
    material::phytoplankton p;
    p.set_wavelength(550e-9);
    p.chl_concentration(10e-6);
    check_close(p.scattering_coefficient(),1,30);
    p.chl_concentration(0.1e-6);
    check_close(p.scattering_coefficient(),0.06,30);
    p.set_wavelength(440e-9);
    p.chl_concentration(5e-6);
    check_close(p.absorption_coefficient(),0.1,30);
    p.chl_concentration(25e-6);
    check_close(p.absorption_coefficient(),0.3,30);
    auto p_forward = p.mueller_matrix(unit_vector{0,0}).value(0,0);
    auto unit_back = unit_vector{constants::pi,0};
    auto p_backward = p.mueller_matrix(unit_back).value(0,0);
    check(p_forward/p_backward > 1e4);
  } end_test_case()
}
