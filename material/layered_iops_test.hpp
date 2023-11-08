#include "layered_iops.hpp"
#include "henyey_greenstein.hpp"
#include "water/cdom.hpp"
#include "mixture.hpp"
#include "atmosphere_ocean.hpp"

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
    check_close(l.alpha_terms(0)[0][0],1/(4*constants::pi));
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
    auto m = std::make_shared<material::atmosphere_ocean>();
    layered_iops iops(m,{-10000,-10,1e-6,0,10,100e3},20);
    m->set_position({0,0,-10});
    iops.set_wavelength(400e-9);
    size_t ly = 0;
    double g = iops.alpha_terms(0)[ly][1]*4*constants::pi/3;
    check_close(g,0.35,50);

  } end_test_case()
}
