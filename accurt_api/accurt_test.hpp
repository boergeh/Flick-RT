#include "../environment/unit_test.hpp"
#include "../material/atmosphere.hpp"
#include "../material/layered_iops.hpp"
#include "../material/henyey_greenstein.hpp"
#include "../material/atmosphere_ocean.hpp"
#include "accurt.hpp"

namespace flick {
  begin_test_case(accurt_test_A) {
    stdvector layer_boundaries = {1, 10e3,20e3,100e3};
    size_t n_terms = 10;
    material::atmosphere::configuration c;
    c.set<size_t>("angles",30);
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
    absorption_coefficient a{0};
    scattering_coefficient b{0.5};
    asymmetry_factor g{0};
 
    accurt::configuration c;
    c.set<std::string>("DETECTOR_ORIENTATION","down");
    c.set<double>("REFERENCE_DETECTOR_HEIGHT",1);
    c.set<double>("BOTTOM_BOUNDARY_SURFACE_SCALING_FACTOR",0);
    auto m = std::make_shared<material::henyey_greenstein>(a,b,g);
    auto ac =  accurt(c,m);
    // van de Hulst 1980, vol 1, chapter 9, table 12, p258, FLUX
    check_close(ac.relative_radiation().y()[0],0.34133, 0.005_pct);
  } end_test_case()

  begin_test_case(accurt_test_D) {
    accurt::configuration ac;
    ac.set<std::string>("DETECTOR_ORIENTATION","up");
    ac.set<double>("REFERENCE_DETECTOR_HEIGHT",100e3);
    ac.set<double>("DETECTOR_HEIGHT",1);
    //material::atmosphere_ocean::configuration mc;
    material::atmosphere::configuration mc;
    mc.set<size_t>("angles",30);
    mc.set<size_t>("heights",8);
    //std::cout << mc;
    //auto m = std::make_shared<material::atmosphere_ocean>(mc);
    auto m = std::make_shared<material::atmosphere>(mc);
    auto a =  accurt(ac,m);
    check_close(a.relative_radiation().y()[0],0.3,50_pct);
    ac.set<double>("DETECTOR_HEIGHT",-0.5);
    auto a2 =  accurt(ac,m);
    check_close(a2.relative_radiation().y()[0],0.3,50_pct);
    //std::cout << a2.relative_radiation() << std::endl;
    ac.set<double>("REFERENCE_DETECTOR_HEIGHT",-0.1);
    auto a3 =  accurt(ac,m);
    check_close(a3.relative_radiation().y()[0],1,10_pct);
    //std::cout << a3.relative_radiation() << std::endl;
  } end_test_case()
}
