#include "atmosphere.hpp"
#include "layered_iops.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    using namespace material;
    atmosphere::configuration c;
    c.set<size_t>("n_angles",50);
    c.set<size_t>("n_heights",6);
    auto atm = std::make_shared<atmosphere>(c);
    auto wls = {400e-9, 800e-9};
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
    auto f = phase_function(*atm);  
  } end_test_case()
}
