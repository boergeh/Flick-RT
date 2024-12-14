#include "../numeric/units.hpp"
#include "ocean.hpp"

namespace flick { 
    begin_test_case(ocean_test_A) {
    using namespace material;
    ocean::configuration c;
    c.set<size_t>("n_angles",100);
    auto oce = ocean{c};
    size_t n1 = oce.material_ids().size();
    c.set<std::string>("mp_names",{"SD16_VF17","SD16_VF18"});
    c.set<double>("mp_concentrations",{1e-3,1e-3});
    c.set<double>("mp_scattering_scaling_factors",{1,1});
    c.set<double>("mp_bleaching_factors",{0,0});
    oce = ocean{c};
    check(oce.material_ids().size()==n1+2);
  } end_test_case()
  
  begin_test_case(ocean_test_B) {
    using namespace material;
    ocean::configuration c;
    double d = 200;
    double wl = 390e-9;
    double distance = d;
    c.set<double>("bottom_depth",d);
    c.set<double>("cdom_440",10);
    auto oce = ocean{c};
    oce.set_wavelength(wl);
    oce.set_position({0,0,0});
    oce.set_direction({0,0,-1});
    double od1 = oce.absorption_optical_depth(distance);
    c.set<double>("concentration_relative_depths",{0,0.5,0.5001,1});
    c.set<double>("concentration_scaling_factors",{0,0,1,1});
    oce = ocean{c};
    oce.set_wavelength(wl);
    oce.set_position({0,0,0});
    oce.set_direction({0,0,-1});
    double od2 = oce.absorption_optical_depth(distance);
    check_close(od1/od2,2,0.1_pct);
  } end_test_case()
  
   begin_test_case(ocean_test_C) {
    // Check no throw for negative abs coefs. at 760 nm
    using namespace material;
    ocean::configuration c;
    c.set<size_t>("n_angles",50);
    c.set<std::string>("mp_names","ECOSENS_HF22_D1");
    c.set<double>("mp_concentrations",1e-3);
    auto oc = std::make_shared<ocean>(c);
    oc->set_wavelength(760e-9);
  } end_test_case()
}
