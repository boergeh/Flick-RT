#include "atmosphere_ocean.hpp"
#include "../numeric/units.hpp"

namespace flick {
  begin_test_case(atmosphere_ocean_test) {
    using namespace units;
    using namespace material;
   
    atmosphere_ocean::configuration c;
    double depth = 100;
    c.set<size_t>("n_angles",50);
    c.set<double>("aerosol_od",1e-5);
    c.set<double>("bottom_depth",depth);
    auto ao = std::make_shared<atmosphere_ocean>(c);
    ao->set_wavelength(400e-9);
    ao->set_position({0,0,0.1});
    double s1 = ao->scattering_coefficient();
    ao->set_wavelength(800e-9);
    double s2 = ao->scattering_coefficient();
    check_close(s2/s1, pow(1./2,4), 10_pct);
    ao->set_position({0,0,-depth});
    ao->set_wavelength(800e-9);
    double a2 = ao->absorption_coefficient();
    check_close(a2, 2, 10_pct);
    ao->set_wavelength(800e-9);
    layered_iops layered_all(ao,{-300,-10,-0.001,0,100e3}, 10);
    layered_all.set_wavelength(800e-9);
    check(layered_all.absorption_coefficient()[0]>1);
    check(layered_all.absorption_coefficient()[1]>1);

  } end_test_case()
}
