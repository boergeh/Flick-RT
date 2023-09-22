#include "../environment/unit_test.hpp"
#include "../material/atmosphere.hpp"
#include "../material/layered_iops.hpp"
#include "../material/henyey_greenstein.hpp"
#include "accurt.hpp"

namespace flick {
  begin_test_case(accurt_test_A) {
    stdvector layer_boundaries = {1, 10e3,100e3};
    size_t n_terms = 16;
    material::atmosphere::config c;
    c.set<size_t>("angles",10);
    c.set<size_t>("heights",8);
    auto atm = material::atmosphere(c);
    layered_iops layered_atmosphere(atm,layer_boundaries,n_terms);
    //std::cout << layered_atmosphere;
  } end_test_case()

   begin_test_case(accurt_test_B) {
     material::white_isotropic m(1.0);
     stdvector boundaries{1,2,10,20};
     size_t n_terms = 4;
     layered_iops iops{m, boundaries, n_terms};
     stdvector wavelengths{300e-9, 500e-9};
     accurt_user_specified accurt{iops, wavelengths};
     //std::cout << accurt;
  } end_test_case()
  
  begin_test_case(accurt_test_C) {
    accurt::configuration ac;
    material::atmosphere::config mc;
    mc.set<size_t>("angles",20);
    mc.set<size_t>("heights",8);
    auto m = std::make_shared<material::atmosphere>(mc);
    auto a =  accurt(ac,m);
    std::cout << std::setprecision(5)<<a.relative_radiation();
  } end_test_case()
}
