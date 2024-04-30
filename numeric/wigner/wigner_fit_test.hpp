#include "wigner_fit.hpp"
#include "../physics_function.hpp"
#include "../constants.hpp"

namespace flick {
  begin_test_case(wigner_fit_test_A) {
    double x = -1;
    size_t n_terms = 33;
    auto f = henyey_greenstein(0.6);
    wigner_fit wf(f,0,0,n_terms,fit::relative);
    check_close(2*wf.coefficients().at(0),1/(2*constants::pi),0.1_pct);
    check_close(wf.value(x),f.value(x),0.1_pct);
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
    check_close(wf.value(x),f().value(x));
    stdvector v = wigner_evaluate(wf.coefficients(),{x},m,n);
    check_close(wf.value(x),v[0]); 
  } end_test_case()
}
