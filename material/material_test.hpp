#include "henyey_greenstein.hpp"
#include "tabulated.hpp"
#include "z_profile.hpp"
#include "../numeric/range.hpp"
#include "aerosols/aerosols.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;

  begin_test_case(material_test_A) {   
    double pi = constants::pi;
    double g = 0.9;
    tabulated_phase_function p{hg_phase_function(g,100)};
    check_close(2*pi*p.integral(),1,0.03_pct);
    check_close(p.asymmetry_factor(),g,0.03_pct);
    
    material::tabulated tab(absorption_coefficient{1},
    			    scattering_coefficient{1},
    			    p);

    const std::string path_{"/material"};
    tabulated_phase_function p2 = read<pe_function>(path_+"/tabulated.txt"); 
    check_close(2*pi*p2.integral(),1,0.3_pct);
  } end_test_case()
  
  begin_test_case(material_test_B) {
    material::rural_aerosols ra;
    material::aggregate_z_profile ag({0,1.5e3,100e3},{0,1,3.14});
    ag.set(pose{{0,0,0},unit_vector{0,0,1}});
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
