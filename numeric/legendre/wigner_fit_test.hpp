#include "wigner_fit.hpp"

namespace flick {
  begin_test_case(wigner_fit_test_A) {
    struct f {
      double g = 0.6;
      double value(double mu) const {
	double theta = acos(mu);
	return henyey_greenstein(g).phase_function(theta);
      }
    };
    size_t n_terms = 16;
    wigner_fit wf(f(),0,0,n_terms,fit::relative);
    check_close(2*wf.coefficients().at(0),1/(2*constants::pi),0.1,"a");
    check_close(wf.value(-1),f().value(-1),0.1,"b");
 
    auto p = flick::read<flick::pl_function>("./petzold_phase_function.txt");
    n_terms = 30;
    wigner_fit wf_p(p,0,0,n_terms,fit::relative);
    check_close(wf_p.value(-1),p.value(-1),0.2,"c");
  } end_test_case()
  
  begin_test_case(wigner_fit_test_B) {
    struct f {
      double value(double x) const {
	return (x+1)*(x-1)*x*x;
      }
    };
    size_t n_terms = 16;
    wigner_fit wf(f(),0,2,n_terms,fit::absolute);
    check_close(wf.value(0.1),f().value(0.1),1e-13);    
  } end_test_case()
}
