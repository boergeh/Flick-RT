#include "layered_iops.hpp"
#include "henyey_greenstein.hpp"
#include "z_profile.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;

  begin_test_case(layered_iops_test) {   
    material::white_isotropic m(1.0);
    m.set(pose{{0,1,1.5},{pi,0}});
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
     m.set(pose{{0,0,10},{pi,0}});
     stdvector bottoms{1,2,10};
     size_t n_terms{4};
     layered_iops iops{m, bottoms, n_terms};
     stdvector wavelengths{300e-9, 500e-9};
     accurt_user_specified accurt{iops, wavelengths};

     stdvector h{0,1,2,10};
     stdvector ang{1,2,3.14};
     material::aggregate a(ang);
     //material::generated_z_profile mz(m,h,ang);
     //material::z_profile mz(m,h);
     //material::aggregate_z_profile ma(mz, h,ang);
     bottoms = {10};
     // m.set(m.pose().move_by({0,0,2}));
     std::cout << accurt;
  } end_test_case()
}
