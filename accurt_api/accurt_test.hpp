#include "../environment/unit_test.hpp"
#include "../material/atmosphere.hpp"
#include "../material/layered_iops.hpp"
#include "accurt.hpp"

namespace flick {
  begin_test_case(accurt_test_A) {
    stdvector layer_boundaries = {1, 10e3,100e3};
    size_t n_terms = 16;
    material::atmosphere::config c;
    c.set<size_t>("angles",100);
    c.set<size_t>("heights",8);
    auto atm = std::make_shared<material::atmosphere>(c);
    layered_iops layered_atmosphere(*atm,layer_boundaries,n_terms);
    //std::cout << layered_atmosphere;
  } end_test_case()
  
  begin_test_case(accurt_test_B) {
    accurt_configuration ac;
    material::atmosphere::config mc;
    mc.set<size_t>("angles",100);
    mc.set<size_t>("heights",8);
    auto m = std::make_shared<material::atmosphere>(mc);
    auto a =  accurt(ac,m);
    //std::cout << a.relative_irradiance();
  } end_test_case()
}
