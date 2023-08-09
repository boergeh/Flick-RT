#include "legendre.hpp"
#include "../physics_function.hpp"
#include "../constants.hpp"

namespace flick {
  struct f {
    double g = 0.3;
    double value(double mu) const {
      return henyey_greenstein(g).value(mu);
    }
  };
  begin_test_case(legendre_test_A) {
    auto p = legendre(129,{0.5});
    check_close(p.value(64,0), -0.0755712307, 1e-8_pct);
    check_close(p.value(128,0), -0.01953466424, 1e-8_pct);
  } end_test_case()
  
  begin_test_case(legendre_test_B) {
    size_t n_terms = 12;
    size_t log2_n_points = 7;
    check_close(gl_integral(std::make_shared<f>(),log2_n_points).value(-1,1),
    		1/(2*constants::pi),0.2_pct);
    std::vector<double> terms = legendre_expansion(std::make_shared<f>(),n_terms,log2_n_points); 
    for (size_t i=0; i < terms.size(); ++i) {
      check_close(terms[i], pow(f().g,i)*(2.*i+1)/2/(2*constants::pi), 0.01_pct);
    }
    std::vector<double> x = read_quadrature(log2_n_points).column(0);
    std::vector<double> v = legendre_evaluation(terms).values(x);
    for (size_t i=0; i < v.size(); ++i) {
      check_close(v[i], f().value(x[i]), 0.02_pct);
    }
  } end_test_case()

   begin_test_case(legendre_test_accumulated_integral) {
    double p = 1e-12;
    accumulated_integral ai(std::make_shared<f>(),p);
    ai.keep_integration_points(true);
    check_close(ai.add_and_get_total(-1,1),1/(2*constants::pi),p);
    //std::cout << std::setprecision(7)<<ai.integration_points();

    p = 0.1;
    double x0 = 1e-6;
    auto d = std::make_shared<log_normal_distribution>(log(x0),0.5);
    accumulated_integral ai2(d,p);
    ai2.keep_integration_points(true);
    for (size_t i=0; i<2; i++) {
      double x1 = x0;
      double x2 = x0;
      ai2.reset_convergence();
      while(not ai2.has_converged_in_one_iteration()) {
	if (i==0) {
	  x2 = x1*0.5;
	  ai2.add_value(x2,x1);
	} else {
	  x2 = x1*2;
	  ai2.add_value(x1,x2);
	}
	x1 = x2;
      }
    }
    check_close(ai2.get_total(),1,p);
    //std::cout << std::setprecision(7)<<ai2.integration_points();
  } end_test_case()
}
