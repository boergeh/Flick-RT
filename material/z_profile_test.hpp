#include "z_profile.hpp"
#include "aerosols/aerosols.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;
 
  begin_test_case(z_profile_aggregate_test) {
    material::rural_aerosols ra;
    material::aggregate ag({0,1,3.14}, {0,1e3,100e3});
    ag.set_direction({0,0,1});
    ag.add(ra);
    check_close(ag.scattering_optical_depth(100e3),
		ra.scattering_optical_depth(100e3));
    unit_vector u{0,0};
    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0));
    ag.add(ra);
    check_close(ag.absorption_optical_depth(100e3),
    		2*ra.absorption_optical_depth(100e3));
    ag.add(ag);
    check_close(ag.absorption_optical_depth(100e3),
    		4*ra.absorption_optical_depth(100e3));

    u = unit_vector{3.14,0};
    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0));
    material::urban_aerosols ua;
    ag.add(ua);
    check_close(ag.absorption_optical_depth(100e3),
    		4*ra.absorption_optical_depth(100e3)+
		ua.absorption_optical_depth(100e3));

    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0)*0.999,0.07_pct); 
  } end_test_case()
}
