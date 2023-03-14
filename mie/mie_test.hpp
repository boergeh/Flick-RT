#include "mie.hpp"

namespace flick {
  begin_test_case(mie_test) {
    parameterized_monodisperesed_mie pmm({1,0},{1.33,0},500e-9);
    double r = 0.5;
    pmm.radius(r);
    pmm.precision(4);
    check_small(pmm.absorption_cross_section(), 1e-12);
    check_close(pmm.extinction_cross_section(), 2*3.14159*pow(r,2), 1);
    check(pmm.scattering_function(0,0)[0] > pmm.scattering_function(0,0)[1]);
 
    polydispersed_mie pm(pmm,log_normal_distribution{log(r),1e-5});    
    check_close(pm.extinction_cross_section(),
    		pmm.extinction_cross_section(),0.01);
    
  } end_test_case()
}
