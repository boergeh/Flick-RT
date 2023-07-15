#include "henyey_greenstein.hpp"
#include "tabulated.hpp"
#include "../environment/input_output.hpp"

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
}
