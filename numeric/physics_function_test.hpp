#include "physics_function.hpp"

namespace flick {
  begin_test_case(physics_function_test_A) {
    using namespace constants;
    henyey_greenstein hg{0};
    double P_iso = 1/(4*pi);
    check_close(hg.phase_function(0), P_iso);
    check_small(hg.inverted_accumulated_angle(0));
    check_close(hg.inverted_accumulated_angle(1), pi);
    hg = henyey_greenstein{0.9};
    check_small(hg.inverted_accumulated_angle(0),1e-7);
    check_close(hg.inverted_accumulated_angle(1), pi,1e-6_pct);
  } end_test_case()
  
  begin_test_case(physics_function_test_B) {
    using namespace constants;
    fournier_forand ff(4.75,1.3);
    std::vector<double> mu = range(-1,1,1000).linspace();
    mu.pop_back();
    std::vector<double> p;
    for (size_t i = 0; i<mu.size(); i++) {
      p.push_back(ff.value(mu[i]));
    }
    pl_function f(mu,p);
    check_close(2*pi*f.integral(-1,1),1,0.2);
  } end_test_case()
  
  begin_test_case(physics_function_test_C) {
    using namespace constants;
    double g0 = 0.5;
    fournier_forand ff(g0);
    std::vector<double> mu = range(-1,1-1e-4,1000).linspace();
    std::vector<double> p;
    for (size_t i = 0; i<mu.size(); i++) {
      p.push_back(ff.value(mu[i])*mu[i]);
    }
    pl_function f(mu,p);
    double g = 2*pi*f.integral(-1,1);
    check_close(g,g0,3_pct);
  } end_test_case()
}
