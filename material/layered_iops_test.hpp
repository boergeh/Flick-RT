#include "layered_iops.hpp"
#include "henyey_greenstein.hpp"
#include "mixture.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;

  begin_test_case(layered_iops_test_A) {
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

   begin_test_case(layered_iops_test_B) {
     material::white_isotropic m(1.0);
     stdvector ang{1,2,3.14};
     stdvector h{0,1,2,10};
     material::mixture a(ang,h);
     
     m.set_position({0,0,10});
     m.set_direction({0,0,-1});
     stdvector boundaries{1,2,10,20};
     size_t n_terms{4};
     layered_iops iops{m, boundaries, n_terms};
     stdvector wavelengths{300e-9, 500e-9};
     accurt_user_specified accurt{iops, wavelengths};
     //std::cout << accurt;
  } end_test_case()
}
