#include "single_layer_slab.hpp"
#include "../material/henyey_greenstein.hpp"

namespace flick {
  
  percentage p{10};

  const double pi = constants::pi;
  begin_test_case(single_layer_slab_test_A) {
    using namespace flick;
    thickness h{1};
    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0.0};
    model::single_layer_slab slab{h};
    slab.fill<material::henyey_greenstein>(a,b,g);    
    slab.adjust_accuracy(percentage{20});
    check_small(slab.hemispherical_reflectance(),1e-12,"a");

    slab.set_bottom(albedo{1});
    check_close(slab.hemispherical_reflectance(),1,1e-12,"b");

    slab.set_bottom(albedo{0.5});
    check_close(slab.hemispherical_reflectance(),0.5,20,"c");
    
    a = 1;
    slab.fill<material::henyey_greenstein>(a,b,g);
    slab.set_bottom(albedo{0});
    check_close(slab.hemispherical_transmittance(),exp(-h()),p(),"d");

  } end_test_case()
  

  begin_test_case(single_layer_slab_test_B) {
    using namespace flick;
    absorption_coefficient a{0};
    scattering_coefficient b{1};
    asymmetry_factor g{0};
    thickness h{1};
    model::single_layer_slab slab{h};
    slab.fill<material::henyey_greenstein>(a,b,g);    
    slab.set_bottom(albedo{0});
    slab.orient_source(zenith_angle{0});
    slab.adjust_accuracy(p);

    // van de Hulst 1980, vol 1, chapter 9, table 12, p258, FLUX
    check_close(slab.hemispherical_reflectance(),0.34133, 1.5*p(),"a");

    // van de Hulst 1980, vol 1, chapter 9, table 12, p259, FLUX
    check_close(slab.hemispherical_transmittance(),0.65867, p(),"b");
  } end_test_case()

  begin_test_case(single_layer_slab_test_C) {
    using namespace flick;
    absorption_coefficient a{0};
    scattering_coefficient b{1};
    asymmetry_factor g{0.5};
    thickness h{1};
    model::single_layer_slab slab{h};
    slab.fill<material::henyey_greenstein>(a,b,g);    
    slab.set_bottom(albedo{0});
    slab.orient_source(zenith_angle{0});
    slab.adjust_accuracy(p);

    // van de Hulst 1980, vol 2, chapter 13, table 35, p419, FLUX
    check_close(slab.hemispherical_transmittance(),0.82389, p());

  } end_test_case()
  
  begin_test_case(single_layer_slab_test_D) {
    using namespace flick;
    double ssalb = 0.9;
    double od = 1;
    double ac = od*(1-ssalb);
    double sc = ssalb*od;
    absorption_coefficient a{ac};
    scattering_coefficient b{sc};
    asymmetry_factor g{0.5};
    thickness h{1};
    model::single_layer_slab slab{h};
    slab.fill<material::henyey_greenstein>(a,b,g);    
    slab.set_bottom(albedo{0});
    slab.orient_source(zenith_angle{0});
    slab.adjust_accuracy(p);
    
    // van de Hulst 1980, vol 2, chapter 13, table 35, p419, FLUX
    check_close(slab.hemispherical_transmittance(),0.73909,p());

  } end_test_case()
  
  begin_test_case(single_layer_slab_test_E) {
    using namespace flick;
    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0.0};
    thickness h{1};
    model::single_layer_slab slab{h};
    slab.fill<material::henyey_greenstein>(a,b,g);    
    slab.set_bottom(albedo{1});
    slab.orient_source(zenith_angle{0});
    slab.adjust_accuracy(p);
    polar_angle pa{0};
    azimuth_angle aa{0};
    vertex_angle acceptance_ang{pi/2};
    double r1 = slab.relative_radiance(pa,aa,acceptance_ang,unit_interval{0});
    check_close(r1,1/pi,p());
  } end_test_case()
  
    begin_test_case(single_layer_slab_test_F) {
    using namespace flick;
    absorption_coefficient a{0};
    scattering_coefficient b{1};
    asymmetry_factor g{0.5};
    thickness h{1};
    model::single_layer_slab slab{h};
    slab.fill<material::henyey_greenstein>(a,b,g);    
    slab.set_bottom(albedo{1});
    slab.orient_source(zenith_angle{0.0});
    slab.adjust_accuracy(percentage{20});
    double r = slab.hemispherical_reflectance(unit_interval{0.5});
    check(r>1,"internal reflectance should be larger than 1 inside white slab");
  } end_test_case()
 
    begin_test_case(single_layer_slab_test_G) {
    const double pi = constants::pi;
    using namespace flick;
    model::single_layer_slab slab{thickness{1}};

    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0.0};
    real_refractive_index n{2.0};

    double brewster = atan(n());

    slab.fill<material::henyey_greenstein>(a,b,g,n());    
    slab.orient_source(zenith_angle{brewster});
    slab.initiate_source(stokes::unpolarized());
    slab.adjust_accuracy(percentage{20});
    slab.set_bottom(albedo{0});
    
    polar_angle theta{brewster};
    azimuth_angle phi{pi};
    vertex_angle acceptance{pi};
    
    double r1 = slab.hemispherical_reflectance();
    check(r1 > 0.01,"a");

    slab.initiate_source(stokes::s_polarized());
    double r2 = slab.hemispherical_reflectance();
    check(r2 > 0.01,"b");
    
    slab.initiate_source(stokes::p_polarized());
    double r3 = slab.hemispherical_reflectance();
    check_small(r3,1e-6,"c");  
  } end_test_case()

   begin_test_case(single_layer_slab_test_H) {
    // Exam UiB, Feb. 2022, PHYS2005, Problem 2e 
    const double pi = constants::pi;
    using namespace flick;
    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0.0};
    real_refractive_index n{1.33};
    double theta0 = pi/4;
    model::single_layer_slab slab{thickness{1}};
    slab.adjust_accuracy(p);    
    slab.initiate_source(stokes::s_polarized());
    slab.orient_source(zenith_angle{theta0});
    slab.fill<material::henyey_greenstein>(a,b,g,n());  
    double T = slab.hemispherical_transmittance();
    double incident_irradiance = 1e-3/1e-4*cos(theta0);
    double bottom_irradiance = incident_irradiance*T;
    double power_hitting_shell = bottom_irradiance*1e-6;
    check_close(power_hitting_shell,6.7012e-6,p());
  } end_test_case()

   begin_test_case(single_layer_slab_test_I) {
    /*
    const double pi = constants::pi;
    using namespace flick;
    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0.0};
    real_refractive_index n{1/1.33};
    double theta0 = pi/4;
    model::single_layer_slab slab{thickness{1}};
    slab.adjust_accuracy(percentage{5});    
    slab.initiate_source(stokes::s_polarized());
    slab.orient_source(zenith_angle{theta0});
    slab.fill<material::henyey_greenstein>(a,b,g,n());  
    check_close(slab.hemispherical_reflectance(),1,0.0001);
    */
  } end_test_case()
}
