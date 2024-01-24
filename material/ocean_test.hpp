#include "ocean.hpp"
#include "../numeric/units.hpp"

namespace flick {
  begin_test_case(ocean_test_A) {
    using namespace material; 
    ocean::configuration c;
    c.set<size_t>("n_angles",100);
    auto oce = ocean(c);
    oce.set_position({0,0,-199});
    oce.set_wavelength(400e-9);
    double a1 = oce.absorption_coefficient();
    check_close(a1,0.09,10_pct);
    oce.set_wavelength(800e-9);
    double a2 = oce.absorption_coefficient();
    check_close(a2,2,10_pct);
  } end_test_case()
  
    begin_test_case(ocean_test_B) {
    using namespace material;
    ocean::configuration c;
    c.set<size_t>("n_angles",100);
    auto oce = ocean{c};
    size_t n1 = oce.material_ids().size();
    c.set<std::string>("mp_names",{"SD16_VF17","SD16_VF18"});
    c.set<double>("mp_concentrations",{1,1});
    oce = ocean{c};
    check(oce.material_ids().size()==n1+1);
  } end_test_case()
}
