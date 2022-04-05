#include "vector.hpp"
namespace flick {
  begin_test_case(vector_test) {
    vector v = {1,0,0};
    vector v2 = v*4 + vector{-1,0,0}*4;
    check_small(norm(v2),1e-15);

    check_small(sin(cast_to_valarray(v2))[0],1e-15);
    
    std::vector<unit_vector> vs{{pi,0}};
    unit_vector uv{0,0,1.01};
    check_small(uv.r()-1,1e-15);

    //double ddd=10_km;

  } end_test_case()
}
