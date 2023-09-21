#include "ab_functions.hpp"
#include "henyey_greenstein.hpp"

namespace flick {
  begin_test_case(ab_functions_test_A) {
    using namespace constants;
    using namespace flick;
    auto [a,b,x] = material::mueller_ab_functions(material::white_isotropic(1),3);
    check_close(a[0][0],1/(4*constants::pi));
  } end_test_case()
  
  begin_test_case(ab_functions_test_B) {
    using namespace constants;
    using namespace flick;
    double g = 0.5;
    auto hg = material::henyey_greenstein(absorption_coefficient{1},
				   scattering_coefficient{1},
				   asymmetry_factor{g});    
    /*
    auto hg = std::make_shared<material::henyey_greenstein>(absorption_coefficient{1},
				   scattering_coefficient{1},
				   asymmetry_factor{g});    
    */
    material::phase_function pf(hg);
    tabulated_phase_function hgpf{hg_phase_function(g,100)};
    check_close(pf.value(0),hgpf.value(constants::pi/2));
    auto [a, b, x] = material::mueller_ab_functions(hg,100);
    check_close(a[0][0],hgpf.value(constants::pi));
  } end_test_case() 
}
