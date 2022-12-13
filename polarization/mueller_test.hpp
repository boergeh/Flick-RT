#include "algorithm.hpp"
#include "isotropic_mueller.hpp"
#include "rayleigh_mueller.hpp"

namespace flick {
  begin_test_case(mueller_test) {
    using namespace constants;
  
    pl_function f;
    size_t n = 150;
    auto ang = range(0,pi,n).linspace();
    for (size_t i=0; i<n; ++i) {
      mueller rm = rayleigh_mueller(ang[i],0.2);
      f.append({ang[i], rm(0).value * sin(ang[i])});
    }
    check_close(2*pi*f.integral(), 1, 0.01,"b");
    
    angular_mueller am(hg_phase_function(0.9,100));
    check_close(am.value(0,0,0),15.12,1);
    check_small(am.value(1,0,0),1e-12);

    angular_mueller am2(hg_phase_function(0.0,100));
    am2.add(3,3,pl_function{{0,3.14159},{1,1}});
    check_close(am2.value(3,3,0),1,1e-9);
    am.add(am2,0.5);
    check_close(am.value(3,3,0),0.5,1e-9);
  } end_test_case()
}

