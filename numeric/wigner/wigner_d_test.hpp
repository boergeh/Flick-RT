#include "wigner_d.hpp"

namespace flick {
  begin_test_case(wigner_d_test) {
    double p = 1e-13;
    double x = -0.25;
    check_close(wigner_d(x,0,0,1).terms()[0],1,p,"d^0_00");
    check_close(wigner_d(x,0,0,2).terms()[1],x,p,"d^1_00");
    check_close(wigner_d(x,0,0,3).terms()[2],0.5*(3*pow(x,2)-1),p,"d^2_00");

    check_small(wigner_d(x,2,2,1).terms()[0],p,"d^0_22");
    check_small(wigner_d(x,2,2,2).terms()[1],p,"d^1_22");
    check_close(wigner_d(x,2,2,3).terms()[2],0.25*pow(1+x,2),p,"d^2_22");
    check_close(wigner_d(x,2,-2,3).terms()[2],0.25*pow(1-x,2),p,"d^2_2,-2");
    check_close(wigner_d(x,0,2,3).terms()[2],
		0.5*pow(3./2,0.5)*(1-pow(x,2)),p,"d^2_02");
    check_close(wigner_d(x,0,2,4).terms()[3],
		pow(15./8,0.5)*(1-pow(x,2))*x,p,"d^3_02");
  } end_test_case()
  
}
