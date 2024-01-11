#include "layered_iops.hpp"
#include "henyey_greenstein.hpp"
#include "water/cdom.hpp"
#include "mixture.hpp"
#include "atmosphere_ocean.hpp"
#include <numbers>

namespace flick {
  begin_test_case(layered_iops_test_A) {
    auto m = std::make_shared<material::white_isotropic>(1.0);
    check_throw(layered_iops(m,{0},4));
    layered_iops l(m,{0,1},4);
    check_close(l.scattering_optical_depth()[0],1);
  } end_test_case()
 
  begin_test_case(layered_iops_test_B) {
    using namespace constants;
    auto m = std::make_shared<material::white_isotropic>(1.0);
    layered_iops l(m,{0,1,9,19},4);
    check_close(l.scattering_optical_depth()[0],1);
    check_close(l.scattering_optical_depth()[1],8);
    check_close(l.scattering_optical_depth()[2],10);
    check_small(l.absorption_optical_depth()[2]);
    check_close(l.single_scattering_albedo()[2],1);
    check_close(l.alpha_terms(0)[0][0],1/(4*std::numbers::pi));
    check_small(l.alpha_terms(0)[0][1]);
    check_close(l.refractive_index()[2],1);   
  } end_test_case()
 
  begin_test_case(layered_iops_test_C) {
    auto m = std::make_shared<material::cdom>();
    layered_iops iops(m,{-1,0,1},10);
    iops.set_wavelength(400e-9);
    double a1 = iops.absorption_coefficient()[0];
    iops.set_wavelength(700e-9);	
    double a2 = iops.absorption_coefficient()[1];
    check(a1/a2 > 2);
  } end_test_case()

  begin_test_case(layered_iops_test_D) {
    double g0 = 0.6;
    auto m = std::make_shared<material::henyey_greenstein>(1,1,g0);
    size_t n_terms = 20;
    layered_iops iops(m,{-10,1e-6,0},n_terms);
    size_t ly = 0;
    double g = iops.alpha_terms(0)[ly][1]*4*std::numbers::pi/3;
    check_close(g0,g,0.1_pct);
  } end_test_case()
  
  begin_test_case(layered_iops_test_E) {
    auto c = material::atmosphere_ocean::configuration();
    c.set<double>("aerosol_od",1);
    c.set<double>("cdom_440",0);
    c.set<double>("chl_concentration",0);
    c.set<double>("nap_concentration",0);
    c.set<double>("bubble_volume_fraction",0);
    auto m = std::make_shared<material::atmosphere_ocean>(c);
    size_t n_terms = 34;
    layered_iops iops(m,{-1000, -100, -10,1, 1e3},n_terms);
    size_t ly = 0;
    double x0 = iops.alpha_terms(0)[ly][0]*4*std::numbers::pi;
    double g = iops.alpha_terms(0)[ly][1]*4*std::numbers::pi/3;
    double x2 = iops.alpha_terms(0)[ly][2];
    double x3 = iops.alpha_terms(0)[ly][3];
    double epsilon = 1e-5;
    check_close(x0,1,0.01_pct);
    check_small(g,epsilon);
    check(x2>epsilon);
    check_small(x3,epsilon);
  } end_test_case()
}
