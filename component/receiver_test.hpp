#include "receiver.hpp"
namespace flick {
  begin_test_case(receiver_test) {
    receiver re;
    re.activate();
    direction_generator dg;
    while (re.received_packages() < 5) {
      radiation_package rp({{0,0,0},{0,0}},stokes{1,0,0,0});
      rp.traveling_direction(dg.isotropic());
      re.receive(rp);
    }
    //std::cout << re;
				    
  } end_test_case()
}
