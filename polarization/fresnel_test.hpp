#include "fresnel.hpp"

namespace flick {
  begin_test_case(fresnel_test_A) {
    std::complex<double> ref{1.33, 0};
    double R = pow((1-ref.real())/(1+ref.real()), 2);
    fresnel fr(ref,1);
    check_close(fr.R(),R,1e-3_pct);    
  } end_test_case()
  
  begin_test_case(fresnel_test_B) {
    std::complex<double> ref{1, 0};
    fresnel f(ref,1);
    check_small(f.reflection_angle());    
  } end_test_case()
}

