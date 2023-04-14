#include "legendre.hpp"
#include "../physics_function.hpp"
#include "../constants.hpp"

namespace flick {
  begin_test_case(legendre_test_A) {
    auto p = legendre(129,{0.5});
    check_close(p.value(64,0),-0.0755712307,1e-8);
    check_close(p.value(128,0),-0.01953466424,1e-8);
  } end_test_case()
  
  begin_test_case(legendre_test_B) {
    struct f {
      double g = 0.3;
      double value(double mu) const {
	double theta = acos(mu);
	return henyey_greenstein(g).phase_function(theta);
      }
    };
    size_t n_terms = 12;
    size_t quad_points = 128;
    check_close(gl_integral(f(),quad_points).value(-1,1),
    		1/(2*constants::pi),0.2);
    std::vector<double> terms = legendre_expansion(f(),n_terms,quad_points); 
    for (size_t i=0; i < terms.size(); ++i) {
      check_close(terms[i],pow(f().g,i)*(2.*i+1)/2/(2*constants::pi),0.01);
    }
    std::vector<double> x = read_quadrature(quad_points).column(0);
    std::vector<double> v = legendre_evaluation(terms).values(x);
    for (size_t i=0; i < v.size(); ++i) {
      check_close(v[i],f().value(x[i]),0.02);
    }
  } end_test_case()
}
