#include "sorted_vector.hpp"

namespace flick {
  begin_test_case(sorted_vector_test) {
    try {
      sorted_vector sv1{{1,-1}};
      check(false);
    } catch (...) {
      check(true);
    }
    
    sorted_vector sv2{{-1,1,2}};
    check(sv2.find(0)==0);
    check(sv2.find(-2)==0);
    check(sv2.find(3)==1);

    sorted_vector sv3;
    for (size_t n=1; n<1e5; ++n)
      sv3.append(exp(1e-4*n));

    using time = std::chrono::steady_clock;
    auto start1 = time::now();
    sv3.find(90);
    auto end1 = time::now();
    std::chrono::duration<double> t1 = end1-start1;

    auto start1b = time::now();
    sv3.set_step_type(step_type::exponential);
    sv3.find(90);
    auto end1b = time::now();
    std::chrono::duration<double> t1b = end1b-start1b;
    check(t1b < t1);
     
    sorted_vector sv4 = sv3;
    sv4.find(2e-4);
    sv4.log_transform();
    auto start2 = time::now();
    sv4.find(log(90));
    auto end2 = time::now();
    std::chrono::duration<double> t2 = end2-start2;
    check(t1 > t2);

    sorted_vector sv5 = sv4;
    sv5.exp_transform();
    check(sv5[100]==sv3[100]);
  } end_test_case()
}
