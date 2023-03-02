#include "delta_fit.hpp"

namespace flick {
  begin_test_case(delta_fit_test) {
    struct f {
      double g = 0.6;
      double value(double mu) const {
	double theta = acos(mu);
	return henyey_greenstein(g).phase_function(theta);
      }
    };
    size_t n_terms = 16;
    delta_fit df(f(),n_terms);
    check_close(2*df.coefficients().at(0),1/(2*constants::pi),0.1);
    check_close(df.function_values({-1}).at(0),f().value(-1),0.1);
    
  } end_test_case()
}
