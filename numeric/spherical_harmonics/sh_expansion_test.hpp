#include "sh_expansion.hpp"

namespace flick {
  begin_test_case(sh_expansion_test) {
    sh_expansion e(1);
    e.set_coefficient_lm(0,0,1);
    //e.set_coefficient_lm(1,-1,1);
    //e.set_coefficient_lm(1,0,1);
    //e.set_coefficient_lm(1,1,1);
    double Y_00 = 1./2*sqrt(1/std::numbers::pi);
    check_close(e.value({0,0}),Y_00);
    //e.write_distribution(10,10,std::cout);
    
  } end_test_case() 
}
