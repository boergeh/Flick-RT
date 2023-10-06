#include "atmosphere_ocean.hpp"
#include "../numeric/units.hpp"

namespace flick {
  begin_test_case(atmosphere_ocean_test) {
    using namespace units;
    using namespace material;
   
    atmosphere_ocean::configuration c;
    c.set<size_t>("angles",100);
    c.set<size_t>("heights",8);
    c.set<double>("aerosol_od",1e-5);
    auto atm = std::make_shared<atmosphere_ocean>(c);
    auto wls = {400e-9, 800e-9};
    atm->set_wavelength(400e-9);
    double s1 = atm->scattering_coefficient();
    atm->set_wavelength(800e-9);
    double s2 = atm->scattering_coefficient();
    check_close(s2/s1, pow(1./2,4), 10_pct);

    atm->set_position({0,0,-199});
    atm->set_wavelength(400e-9);
    double a1 = atm->absorption_coefficient();
    check_close(a1,0.02,10_pct);
    atm->set_wavelength(800e-9);
    double a2 = atm->absorption_coefficient();
    check_close(a2,2,10_pct);

    atm->set_wavelength(800e-9);
    layered_iops layered_all(atm,{-300,-0.001,0,100e3}, 3);
    check(layered_all.absorption_coefficient()[0]>1);

  } end_test_case()
}
