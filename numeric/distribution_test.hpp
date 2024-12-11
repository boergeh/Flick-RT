#include "distribution.hpp"

namespace flick {
  begin_test_case(distribution_test_A) {
    using namespace distribution;
    double p = 1e-12;
    erf_inv erf_inv(p);
    check_small(erf_inv(0), p);
    double x = 0.99;
    check_close(std::erf(erf_inv(x)), x, p);
    check_close(erf_inv(0.999), 2.326753765513524, p);
    check_close(erf_inv(-0.999), -2.326753765513524, p);
  } end_test_case()
  
  begin_test_case(distribution_test_B) {
    using namespace distribution;
    double a = 8;
    double b = 12;
    double c = 10;
    triangular t(a,b,c);
    check_small(t.pdf(a));
    check_small(t.pdf(b));
    check_small(t.cdf(a));
    check_close(t.cdf(b),1);
    check_close(t.quantile(0),a);
    check_close(t.quantile(1),b);
    check_close(t.quantile(0.5),a+(b-a)/2);
  } end_test_case()
  
  begin_test_case(distribution_test_C) {
    using namespace distribution;
    std::vector<double> x = {-3, -2, -1, 0, 1, 2};
    std::vector<double> y = {9, 4, 1, 0.001, 1, 4};
    pe_function f{x,y};
    //curvature c(f);
    //check_small(curvature(f).cdf(-5));
    //check_close(curvature(f).cdf(5),1);
    //check_close(curvature(f).quantile(0),-3);
  } end_test_case()
}
