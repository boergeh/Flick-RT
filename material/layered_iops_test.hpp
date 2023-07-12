#include "layered_iops.hpp"
#include "henyey_greenstein.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;

  begin_test_case(layered_iops_test) {   
    material::white_isotropic m(1.0);
    m.position({0,1,1.5});
    m.direction({pi,0});
    layered_iops l(m,{1,2,10},4);
    check_close(l.scattering_optical_depth()[0],1);
    check_close(l.scattering_optical_depth()[1],1);
    check_close(l.scattering_optical_depth()[2],8);
    check_small(l.absorption_optical_depth()[2]);
    check_close(l.single_scattering_albedo()[2],1);
    check_close(l.alpha_terms(0)[0][0],1/(4*constants::pi));
    check_small(l.alpha_terms(0)[0][1]);
    check_close(l.refractive_index()[2],1);   
  } end_test_case()

   begin_test_case(accurt_user_specified_test) {
     material::white_isotropic m(1.0);
     layered_iops l(m,{1,2,10},4);
     accurt_user_specified u(l,{300e-9, 500e-9});
     //std::cout << u;
  } end_test_case()
}
