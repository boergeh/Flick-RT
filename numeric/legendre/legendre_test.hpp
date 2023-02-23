#include "legendre.hpp"
#include "../physics_function.hpp"
#include "../constants.hpp"

namespace flick {
  begin_test_case(legendre_test) {
    struct f {
      double g = 0.3;
      double value(double mu) const {
	double theta = acos(mu);
	return henyey_greenstein(g).phase_function(theta);
      }
    };
    size_t n_terms = 12;
    size_t quad_points = 100;
    check_close(gl_integral(f(),quad_points).value(-1,1),
		1/(2*constants::pi),0.2);
    std::vector<double> terms = legendre_expansion(f(),n_terms,quad_points); 
    for (size_t i=0; i < terms.size(); ++i) {
      //std::cout << terms[i] << " ";
      check_close(terms[i],pow(f().g,i)*(2.*i+1)/2/(2*constants::pi),0.01);
    }
    std::vector<double> x = read_quadrature(quad_points).column(0);
    std::vector<double> v = legendre_evaluation(terms).values(x);
    //pl_function hg;
    //f fhg;
    for (size_t i=0; i < v.size(); ++i) {
      //hg.append({x[i],f().value(x[i])});
      check_close(v[i],f().value(x[i]),0.02);
    }
    //write(hg,"/numeric/gl_integral_hg.txt");
    //n_terms = 8;
    //terms = delta_log(f(),n_terms,quad_points);
    //std::cout << std::endl;
    //for (size_t i=0; i < terms.size(); ++i) {
    //  std::cout << terms[i] << " ";
    //}
    //pl_function p{x,legendre_evaluation(terms).values(x)};
    //write(p,"/numeric/gl_integral_delta_log.txt");
    
  } end_test_case()
}
