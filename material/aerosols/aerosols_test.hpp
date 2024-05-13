#include "../../environment/unit_test.hpp"
#include "aerosols.hpp"

namespace flick {
  using namespace flick;
  begin_test_case(aerosols_test_A) {
    material::rural_aerosols ra(300,0.1,0.9);
    ra.set(pose{{0,0,1},unit_vector{0,0,1}});
    ra.set_wavelength(500e-9);
    check_close(ra.absorption_coefficient(),1e-5,40.0_pct);
  } end_test_case()

  begin_test_case(aerosols_test_B) {
    material::rural_aerosols ra(300,0.1,0.5);
    ra.set_wavelength(500e-9);
    size_t n_terms = 15;
    auto [alpha, beta] = material::fitted_mueller_alpha_beta(ra,n_terms);
    double c0 = alpha[0][0];
    double g = alpha[0][1]*4*constants::pi/3;
    check_close(c0,1/(4*constants::pi),10_pct);
    check_close(g,0.6575,0.5_pct);
  } end_test_case()

  begin_test_case(aerosols_test_C) {
    material::rural_aerosols ra(300,0.1,0.5);
    ra.set_wavelength(400e-9);
    double b_400 = ra.scattering_coefficient();
    ra.set_wavelength(700e-9);
    double b_700 = ra.scattering_coefficient();
    check_close(b_400/b_700, 2, 10_pct);
  } end_test_case()
  
  begin_test_case(aerosols_test_D) {
    // Check that scattering slope change with humidity
    material::rural_aerosols low_rh(300,0.1,0);
    low_rh.set_wavelength(400e-9);
    material::rural_aerosols high_rh(300,0.1,1);
    high_rh.set_wavelength(400e-9);
    double r = low_rh.scattering_coefficient()/high_rh.scattering_coefficient();
    check(r > 1.01);
  } end_test_case()
}
