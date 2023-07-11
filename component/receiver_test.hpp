#include "receiver.hpp"

namespace flick {
  begin_test_case(receiver_test) {
    receiver re;
    re.activate();
    direction_generator dg;
    while (re.received_packages() < 5) {
      radiation_package rp({{0,0,0},{0,0}},stokes{1,0,0,0});
      rp.rotate_to(rotation_to(dg.isotropic()));
      re.receive(rp);
    }
  } end_test_case()
}
