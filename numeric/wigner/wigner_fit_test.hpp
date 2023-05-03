#include "wigner_fit.hpp"
#include "../physics_function.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
  begin_test_case(wigner_fit_test_A) {
    struct f {
      double g = 0.6;
      double value(double mu) const {
	double theta = acos(mu);
	return henyey_greenstein(g).phase_function(theta);
      }
    };
    double x = -1;
    size_t n_terms = 16;
    wigner_fit wf(f(),0,0,n_terms,fit::relative);
    check_close(2*wf.coefficients().at(0),1/(2*constants::pi),0.1,"a");
    check_close(wf.value(x),f().value(x),0.1,"b");
 
    auto p = flick::read<flick::pl_function>("./petzold_phase_function.txt");
    n_terms = 49;
    wigner_fit wf_p(p,0,0,n_terms,fit::relative);
    check_close(wf_p.value(x),p.value(x),0.2,"c");
  } end_test_case()
  
  begin_test_case(wigner_fit_test_B) {
    struct f {
      double value(double x) const {
	return (x+1)*(x-1)*pow(x,3);
      }
    };
    size_t n_terms = 6;
    int m = 0;
    int n = 2;
    wigner_fit wf(f(),m,n,n_terms,fit::absolute);
    double x = -0.5;
    double p = 1e-13;
    check_close(wf.value(x),f().value(x),p);
    stdvector v = wigner_evaluate(wf.coefficients(),{x},m,n);
    check_close(wf.value(x),v[0],p);
    
  } end_test_case()
}
