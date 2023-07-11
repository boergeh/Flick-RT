#include "wigner_d.hpp"

namespace flick {
  begin_test_case(wigner_d_test) {
    double x = -0.25;
    check_close(wigner_d(x,0,0,1).terms()[0],1);
    check_close(wigner_d(x,0,0,2).terms()[1],x);   
    check_close(wigner_d(x,0,0,3).terms()[2],0.5*(3*pow(x,2)-1));
    check_small(wigner_d(x,2,2,1).terms()[0]);
    check_small(wigner_d(x,2,2,2).terms()[1]);
    check_close(wigner_d(x,2,2,3).terms()[2],0.25*pow(1+x,2));
    check_close(wigner_d(x,2,-2,3).terms()[2],0.25*pow(1-x,2));
    check_close(wigner_d(x,0,2,3).terms()[2],0.5*pow(3./2,0.5)*(1-pow(x,2)));
    check_close(wigner_d(x,0,2,4).terms()[3],pow(15./8,0.5)*(1-pow(x,2))*x);
  } end_test_case() 
}
