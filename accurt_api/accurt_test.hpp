#include "../environment/unit_test.hpp"
#include "../material/atmosphere.hpp"
#include "../material/layered_iops.hpp"
#include "../material/henyey_greenstein.hpp"
#include "../material/atmosphere_ocean.hpp"
#include "accurt.hpp"

namespace flick {
  begin_test_case(accurt_test_A) {
    // Compare with van de Hulst 1980, vol 1, chapter 9, table 12,
    // p258, FLUX
    absorption_coefficient a{0};
    scattering_coefficient b{0.5};
    asymmetry_factor g{0};
    accurt::configuration c;
    c.set<std::string>("detector_orientation","down");
    c.set<std::string>("detector_type","plane_irradiance");
    c.set<double>("detector_height",1);
    c.set<double>("reference_detector_height",1);
    c.set<size_t>("stream_upper_slab_size",8);
    c.set<double>("bottom_boundary_surface_scaling_factor",0);
    c.set<double>("detector_wavelengths",400e-9);
    auto m = std::make_shared<material::henyey_greenstein>(a,b,g);
    auto ac =  accurt(c,m);
    check_close(ac.relative_radiation().y()[0],0.34133, 0.005_pct);
  } end_test_case()
  
  begin_test_case(accurt_test_B) {
    // Compare with van de Hulst 1980, vol 2, chapter 13, table 35,
    // p431, FLUX
    double ssalb = 0.9;
    double od = 1;
    double ac = od*(1-ssalb);
    double sc = ssalb*od;
    absorption_coefficient a{ac/2}; // 2 m slab thickness
    scattering_coefficient b{sc/2};
    asymmetry_factor g{0.75};
    size_t n_angles = 150;
    accurt::configuration c;
    c.set<std::string>("detector_orientation","up");
    c.set<std::string>("detector_type","plane_irradiance");
    c.set<double>("detector_height",-1);
    c.set<double>("reference_detector_height",1);
    c.set<size_t>("stream_upper_slab_size", c.to_streams(n_angles));
    c.set<double>("bottom_boundary_surface_scaling_factor",0);
    c.set<double>("detector_wavelengths",400e-9);
    auto m = std::make_shared<material::henyey_greenstein>(a,b,g);
    auto acc =  accurt(c,m);
    check_close(acc.relative_radiation().y()[0],0.82792, 0.001_pct);
  } end_test_case()
  
  begin_test_case(accurt_test_C) {
    // Assert low remote sensing reflectance in NIR and that detectors
    // can have same height
    size_t n_angles = 30;
    accurt::configuration ac;
    ac.set<size_t>("stream_upper_slab_size",ac.to_streams(n_angles));
    ac.set<double>("detector_wavelengths",950e-9);
    ac.set<std::string>("detector_orientation","down");
    ac.set<std::string>("detector_type","radiance");
    ac.set<double>("detector_height",0.1);
    ac.set<double>("reference_detector_height",0.1);
    ac.set<std::string>("subtract_specular_radiance","true");
    material::atmosphere_ocean::configuration mc;
    mc.set<size_t>("n_angles",n_angles);
    mc.set<size_t>("n_heights",3);
    mc.set<double>("aerosol_od",1);
    auto m = std::make_shared<material::atmosphere_ocean>(mc);
    auto a =  accurt(ac, m);
    double Rrs = a.relative_radiation().y()[0];
    check_small(Rrs, 0.0002);
  } end_test_case()
  
  begin_test_case(accurt_test_D) {
    // Check increase in nadir radiance below surface
    size_t n_angles = 30;
    accurt::configuration ac;
    ac.set<size_t>("stream_upper_slab_size",ac.to_streams(n_angles));
    ac.set<double>("detector_wavelengths",500e-9);
    ac.set<std::string>("detector_orientation","down");
    ac.set<std::string>("detector_type","radiance");
    ac.set<double>("detector_height",0.01);
    ac.set<double>("reference_detector_height",120e3);
    material::atmosphere_ocean::configuration mc;
    mc.set<size_t>("n_angles",n_angles);
    mc.set<size_t>("n_heights",3);
    auto m = std::make_shared<material::atmosphere_ocean>(mc);
    auto a_above =  accurt(ac,m);
    double L_above = a_above.relative_radiation().y()[0];
    ac.set<double>("detector_height",-0.01);
    auto a_below =  accurt(ac,m);
    double L_below = a_below.relative_radiation().y()[0];
    check_close(L_above, L_below/pow(1.33,2), 3_pct);
  } end_test_case()
  
  begin_test_case(accurt_test_E) {
    // Assert flick default atmosphere-ocean remote sensing reflectance
    size_t n_angles = 30;
    accurt::configuration ac;
    ac.set<size_t>("stream_upper_slab_size",ac.to_streams(n_angles));
    ac.set<double>("detector_wavelengths",500e-9);
    ac.set<std::string>("detector_orientation","down");
    ac.set<std::string>("detector_type","radiance");
    ac.set<double>("detector_height",0.01);
    ac.set<double>("reference_detector_height",120e3);
    ac.set<std::string>("subtract_specular_radiance","true");
    material::atmosphere_ocean::configuration mc;
    mc.set<size_t>("n_angles",n_angles);
    mc.set<size_t>("n_heights",3);
    auto m = std::make_shared<material::atmosphere_ocean>(mc);
    auto a =  accurt(ac, m);
    double Rrs = a.relative_radiation().y().at(0);
    check_close(Rrs, 0.0392, 0.2_pct);
  } end_test_case()
  
  begin_test_case(accurt_test_F) {
    // Check that azimuthally averaged radiance is the same as
    // directional specific radiance at polar angle is close to 180
    size_t n_angles = 200;
    accurt::configuration ac;
    ac.set<double>("source_zenith_angle",60);
    ac.set<size_t>("stream_upper_slab_size",ac.to_streams(n_angles));
    ac.set<double>("detector_wavelengths",500e-9);
    ac.set<std::string>("detector_orientation","down");
    ac.set<std::string>("reference_detector_orientation","up");
    ac.set<std::string>("detector_type","radiance");
    ac.set<double>("detector_height",120e3);
    ac.set<double>("reference_detector_height",120e3);
    material::atmosphere_ocean::configuration mc;
    mc.set<size_t>("n_angles",n_angles);
    mc.set<size_t>("n_heights",3);
    mc.set<double>("aerosol_od",0);
    auto m = std::make_shared<material::atmosphere_ocean>(mc); 
    auto a_avg =  accurt(ac, m);
    double r_avg = a_avg.relative_radiation().y()[0];
    ac.set<double>("detector_orientation_override",{179.91,180});
    auto a_dir =  accurt(ac, m);
    double r_dir = a_dir.relative_radiation().y()[0];
    check_close(r_avg, r_dir, 1.7_pct);
  } end_test_case()
}
