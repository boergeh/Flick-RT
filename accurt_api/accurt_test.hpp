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
    auto atm = std::make_shared<material::atmosphere>(c);
    layered_iops layered_atmosphere(atm,layer_boundaries,n_terms);
    //std::cout << layered_atmosphere;
  } end_test_case()

   begin_test_case(accurt_test_B) {
    auto m = std::make_shared<material::white_isotropic>(1.0);
     stdvector boundaries{1,2,10,20};
     size_t n_terms = 4;
     auto iops = std::make_shared<layered_iops>(m, boundaries, n_terms);
     stdvector wavelengths{300e-9, 500e-9};
     accurt_user_specified accurt{iops, wavelengths};
     //std::cout << accurt;
  } end_test_case()
  
  begin_test_case(accurt_test_C) {
    absorption_coefficient a{0};
    scattering_coefficient b{0.5};
    asymmetry_factor g{0};
 
    accurt::configuration c;
    c.set<std::string>("detector_orientation","down");
    c.set<double>("reference_detector_height",1);
    c.set<double>("BOTTOM_BOUNDARY_SURFACE_SCALING_FACTOR",0);
    auto m = std::make_shared<material::henyey_greenstein>(a,b,g);
    auto ac =  accurt(c,m);
    // van de Hulst 1980, vol 1, chapter 9, table 12, p258, FLUX
    check_close(ac.relative_radiation().y()[0],0.34133, 0.005_pct);
  } end_test_case()
 
  begin_test_case(accurt_test_D) {  
    accurt::configuration ac;
    ac.set<double>("DETECTOR_WAVELENGTHS",400e-9);
    ac.set<std::string>("detector_orientation","up");
    ac.set<double>("reference_detector_height",100e3);
    ac.set<double>("detector_height",1);
    material::atmosphere_ocean::configuration mc;
    mc.set<size_t>("angles",30);
    mc.set<size_t>("heights",8);
    //mc.set<double>("aerosol_od",0);
    //mc.set<double>("cloud_liquid",0);
    //mc.set<double>("pressure",10);

    auto m = std::make_shared<material::atmosphere_ocean>(mc);

    auto a =  accurt(ac,m);
    check_close(a.relative_radiation().y()[0],0.3,50_pct);
    ac.set<double>("detector_height",-0.5);
    auto a2 =  accurt(ac,m);
    check_close(a2.relative_radiation().y()[0],0.3,50_pct);
    ac.set<double>("reference_detector_height",-0.1);
    auto a3 =  accurt(ac,m);
    check_close(a3.relative_radiation().y()[0],1,10_pct);
    ac.set<double>("detector_height",1e-9);
    ac.set<double>("reference_detector_height",-0.9);
    auto a4 =  accurt(ac,m);
    check_close(a4.relative_radiation().y()[0],1,10_pct);
    
    ac.set<double>("DETECTOR_WAVELENGTHS",{400e-9,700e-9});
    ac.set<double>("BOTTOM_BOUNDARY_SURFACE_SCALING_FACTOR",0);
    ac.set<double>("reference_detector_height",120e3);
    ac.set<double>("detector_height",120e3);
    ac.set<std::string>("detector_orientation","down");
    auto a5 =  accurt(ac,m);
    pp_function r = a5.relative_radiation();
    std::cout <<std::setprecision(9)<< r;
    double dr = r.y()[0] / r.y()[1];
    check(dr > 1.01);
  } end_test_case()
  
  begin_test_case(accurt_test_E) {
    configuration_template::toa_reflectance c;
    std::string f = "./toa_reflectance_configuration";
    c.write(f);
    auto c2 = read<configuration_template::toa_reflectance>(f);
    //std::cout << c2;
  } end_test_case()
}
