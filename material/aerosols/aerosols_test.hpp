#include "../../environment/unit_test.hpp"
#include "aerosols.hpp"

namespace flick {
  using namespace flick;
  begin_test_case(aerosols_test_A) {
    material::rural_aerosols ra;
    ra.set(pose{{0,0,1},unit_vector{0,0,1}});
    ra.set_wavelength(500e-9);
    check_close(ra.absorption_coefficient(),1e-5,40.0_pct);
  } end_test_case()

  begin_test_case(aerosols_test_B) {
    material::rural_aerosols ra;
    ra.set_wavelength(500e-9);
    size_t n_terms = 15;
    auto [alpha, beta] = material::fitted_mueller_alpha_beta(ra,n_terms);
    double c0 = alpha[0][0];
    double g = alpha[0][1]*4*constants::pi/3;
    check_close(c0,1/(4*constants::pi),10_pct);
    check_close(g,0.6575,0.5_pct);
  } end_test_case()

}
