#include "../environment/unit_test.hpp"
#include "../material/atmosphere.hpp"
#include "../material/layered_iops.hpp"
#include "../material/henyey_greenstein.hpp"
#include "../material/atmosphere_ocean.hpp"
#include "accurt.hpp"

namespace flick {
  begin_test_case(accurt_test_A) {
    // Compare with van de Hulst 1980, vol 1, chapter 9, table 12,
    // p258, FLUX U and AVERAGE N for mu_0 = 1.0
    absorption_coefficient a{0};
    scattering_coefficient b{0.5};
    asymmetry_factor g{0};
    accurt::configuration c;
    c.set<std::string>("detector_orientation","down");
    c.set<std::string>("detector_type","plane_irradiance");
    c.set<double>("detector_height",1);
    c.set<double>("reference_detector_height",1);
    c.set<size_t>("stream_upper_slab_size",12);
    c.set<double>("bottom_boundary_surface_scaling_factor",0);
    c.set<double>("detector_wavelengths",400e-9);
    auto m = std::make_shared<material::henyey_greenstein>(a,b,g);
    auto ac =  accurt(c,m);
    double FLUX_U = 0.34133;
    check_close(ac.relative_radiation().y()[0],FLUX_U, 0.01_pct);
    c.set<std::string>("detector_type","scalar_irradiance");
    auto ac2 =  accurt(c,m);
    double AVERAGE_N = 0.37869;
    check_close(ac2.relative_radiation().y()[0]/2,AVERAGE_N, 0.01_pct);
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
    // Assert low remote sensing reflectance in red and that detectors
    // can have same height
    size_t n_angles = 40;
    double wl = 620e-9;
    accurt::configuration ac;
    ac.set<size_t>("stream_upper_slab_size",ac.to_streams(n_angles));
    ac.set<double>("detector_wavelengths",wl);
    ac.set<std::string>("detector_orientation","down");
    ac.set<std::string>("detector_type","radiance");
    ac.set<double>("detector_height",0.1);
    ac.set<double>("reference_detector_height",0.1);
    ac.set<std::string>("subtract_specular_radiance","true");
    material::atmosphere_ocean::configuration mc;
    mc.set<size_t>("n_angles",n_angles);
    mc.set<size_t>("n_heights",3);
    mc.set<double>("aerosol_od",0);
    mc.set<double>("water_salinity",35);
    mc.set<double>("water_temperature",273+15);
    auto m = std::make_shared<material::atmosphere_ocean>(mc);
    m->set_wavelength(wl);
    auto a =  accurt(ac, m);
    double Rrs = a.relative_radiation().y()[0];
    check_close(Rrs, 0.11e-3, 5_pct);
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
    mc.set<double>("nap_concentration",1e-3);
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
    size_t n_angles = 40;
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
    mc.set<double>("nap_concentration",1e-3);
    auto m = std::make_shared<material::atmosphere_ocean>(mc);
    auto a =  accurt(ac, m);
    double Rrs = a.relative_radiation().y().at(0);
    check_close(Rrs, 0.0397, 0.3_pct);
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
    mc.set<double>("nap_concentration",1e-3);
    auto m = std::make_shared<material::atmosphere_ocean>(mc); 
    auto a_avg =  accurt(ac, m);
    double r_avg = a_avg.relative_radiation().y()[0];
    ac.set<double>("detector_orientation_override",{179.91,180});
    auto a_dir =  accurt(ac, m);
    double r_dir = a_dir.relative_radiation().y()[0];
    check_close(r_avg, r_dir, 1.7_pct);
  } end_test_case()

   begin_test_case(accurt_test_G) {
    // Check that nadir toa radiance is increasing slightly when adding a
    // thin cloud layer
    size_t n_angles = 100;
    accurt::configuration ac;
    ac.set<size_t>("stream_upper_slab_size",ac.to_streams(n_angles));
    ac.set<double>("detector_wavelengths",400e-9);
    ac.set<std::string>("detector_orientation","down");
    ac.set<std::string>("detector_type","radiance");
    ac.set<double>("detector_height",120e3);
    ac.set<double>("reference_detector_height",120e3);
    ac.set<double>("bottom_boundary_surface_scaling_factor",0);
    material::atmosphere::configuration mc;
    mc.set<size_t>("n_angles",n_angles);
    mc.set<size_t>("n_heights",8);
    mc.set<std::string>("gases","no2");
    auto m_clear = std::make_shared<material::atmosphere>(mc);
    auto a_clear =  accurt(ac,m_clear);
    double L_toa_clear = a_clear.relative_radiation().y()[0];
    mc.set<double>("cloud_liquid",1e-7);
    auto m_cloudy = std::make_shared<material::atmosphere>(mc);
    auto a_cloudy =  accurt(ac,m_cloudy);
    double L_toa_cloudy = a_cloudy.relative_radiation().y()[0];
    check(L_toa_cloudy > L_toa_clear);
    check_close(L_toa_cloudy,L_toa_clear,1_pct);
  } end_test_case()
}
