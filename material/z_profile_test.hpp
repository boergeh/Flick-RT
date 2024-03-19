#include "z_profile.hpp"
#include "henyey_greenstein.hpp"

namespace flick {
  begin_test_case(z_profile_test) {
    auto hg = std::make_shared<material::henyey_greenstein>(1,0,0);
    material::scaled_z_profile<pl_function> p(hg, {-200,-1,0},{0, 0.5, 1});
    p.set_position({0,0,-1});
    check_close(p.absorption_coefficient(),0.5);
    p.set_position({0,0,-200});
    check_small(p.absorption_coefficient());
    check_close(p.mueller_matrix({0,0,1}).value(0,0),1/(4*constants::pi));
  } end_test_case()
}
