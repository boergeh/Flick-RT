#include "atmosphere.hpp"
#include "layered_iops.hpp"

namespace flick {
  begin_test_case(atmosphere_test_A) {
    using namespace material;
    atmosphere::configuration c;
    c.set<size_t>("n_angles",50);
    c.set<size_t>("n_heights",6);
    auto atm = std::make_shared<atmosphere>(c);
    atm->set_wavelength(400e-9);
    double s_a = atm->scattering_coefficient();
    atm->set_wavelength(800e-9);
    double s_b = atm->scattering_coefficient();
    check_close(s_b/s_a, pow(1./2,4),10_pct);
    atm->set_position({0,0,0});
    double s1 = atm->scattering_coefficient();
    auto iops = layered_iops(atm,range(0.1,100e3,8).logspace(),4);
    iops.set_wavelength(500e-9);
    stdvector od = iops.scattering_optical_depth();
    atm->set_position({0,0,0.1});
    check_close(vec::sum(od), atm->scattering_optical_depth(100e3), 0.001_pct);
  } end_test_case()
  
  begin_test_case(atmosphere_test_B) {
    using namespace material;
    double aero_od = 0.1;
    double h_bottom = 0;
    double h_top = 100e3;
    size_t n_terms = 4;
    double wl = 550e-9;
    atmosphere::configuration c;
    auto aero = std::make_shared<rural_aerosols>(h_top,aero_od,0.5);
    stdvector h = range(h_bottom,h_top,3).linspace();
    auto iops = layered_iops(aero, h, n_terms);
    iops.set_wavelength(wl);
    stdvector scat_od = iops.scattering_optical_depth();
    stdvector abs_od = iops.absorption_optical_depth();
    aero->set_position({0,0,h_bottom});
    check_close(vec::sum(scat_od), aero->scattering_optical_depth(h_top), 0.001_pct);
    check_close(vec::sum(abs_od), aero->absorption_optical_depth(h_top), 0.001_pct);
    check_close(vec::sum(abs_od)+vec::sum(scat_od), aero_od, 0.001_pct);
  } end_test_case()
  
  begin_test_case(atmosphere_test_C) {
    using namespace material;
    double aero_od = 0.1;
    double h_bottom = 0.1;
    double h_top = 120e3;
    size_t n_layers = 3;
    size_t n_terms = 4;
    double wl = 550e-9;
    atmosphere::configuration c;
    c.set<size_t>("n_angles",20);
    c.set<size_t>("n_heights",6);
    c.set<double>("pressure",1);
    c.set<double>("aerosol_od",aero_od);
    c.set<double>("ozone",0);
    c.set<double>("water_vapor",0);
    auto atm = std::make_shared<atmosphere>(c);
    material::z_profile<pe_function>* zp = dynamic_cast<material::z_profile<pe_function>*>(&*atm);
    stdvector h = zp->height_grid();
    auto iops = layered_iops(atm, h, n_terms);
    iops.set_wavelength(wl);
    stdvector scat_od = iops.scattering_optical_depth();
    stdvector abs_od = iops.absorption_optical_depth();
    atm->set_position({0,0,h_bottom});
    check_close(vec::sum(scat_od), atm->scattering_optical_depth(h_top), 0.01_pct);
    check_close(vec::sum(abs_od), atm->absorption_optical_depth(h_top), 0.01_pct);
    check_close(vec::sum(abs_od)+vec::sum(scat_od), aero_od, 0.01_pct);
  } end_test_case()
}
