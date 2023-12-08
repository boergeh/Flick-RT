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
    check_close(2*df.coefficients().at(0),1/(2*constants::pi),0.1_pct);
    check_close(df.function_values({-1}).at(0),f().value(-1),0.1_pct);

    auto g = flick::read<flick::pl_function>("./petzold_phase_function.txt");
    n_terms = 30;
    delta_fit df2(g,n_terms);
    check_close(df2.function_values({-1}).at(0),g.value(-1),0.5_pct);
  } end_test_case()
}
