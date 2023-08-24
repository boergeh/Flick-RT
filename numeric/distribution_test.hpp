#include "distribution.hpp"

namespace flick {
  begin_test_case(distribution_test) {
    using namespace distribution;
    double p = 1e-12;
    erf_inv erf_inv(p);
    check_small(erf_inv(0), p);
    double x = 0.99;
    check_close(std::erf(erf_inv(x)), x, p);
    check_close(erf_inv(0.999), 2.326753765513524, p);
    check_close(erf_inv(-0.999), -2.326753765513524, p);

    //std::cout << std::setprecision(5)<<quantiles(log_normal(0,0.1),500);
  } end_test_case()
}
