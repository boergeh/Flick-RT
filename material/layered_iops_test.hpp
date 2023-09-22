#include "layered_iops.hpp"
#include "henyey_greenstein.hpp"
#include "mixture.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;

  begin_test_case(layered_iops_test_A) {
    using namespace constants;
    material::white_isotropic m(1.0);
    check_throw(layered_iops(m,{0},4));
    layered_iops l(m,{0,1},4);
    check_close(l.scattering_optical_depth()[0],1);
    } end_test_case()
 
  begin_test_case(layered_iops_test_B) {
    using namespace constants;
    material::white_isotropic m(1.0);
    layered_iops l(m,{0,1,9,19},4);
    check_close(l.scattering_optical_depth()[0],1);
    check_close(l.scattering_optical_depth()[1],8);
    check_close(l.scattering_optical_depth()[2],10);
    check_small(l.absorption_optical_depth()[2]);
    check_close(l.single_scattering_albedo()[2],1);
    check_close(l.alpha_terms(0)[0][0],1/(4*constants::pi));
    check_small(l.alpha_terms(0)[0][1]);
    check_close(l.refractive_index()[2],1);   
  } end_test_case()
}
