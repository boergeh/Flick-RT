#include "atmosphere.hpp"
#include "../numeric/units.hpp"
#include "water/pure_water.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    using namespace units;
    using namespace material;
   
    atmosphere::config c;
    //std::cout << "heights: "<< c.get<size_t>("heights");
    atmosphere atm{c};
    auto wls = range(280e-9,950e-9,20).logspace();
    auto pw = pure_water();
    //std::cout << std::setprecision(5)
    //	      << optical_depth(atm,100e3,wls).scattering() << std::endl;


  } end_test_case()
}
