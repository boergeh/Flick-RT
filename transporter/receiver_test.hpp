#include "receiver.hpp"
namespace flick {
  begin_test_case(receiver_test) {
    receiver re;
    direction_generator dg;
    while (re.received_packages() < 5) {
      radiation_package rp(600e-9,stokes{1,0,0,0});
      rp.rotate_to(dg.isotropic());
      re.receive(rp);
    }
    //std::cout << re;
				    
  } end_test_case()
}
