#include "ocean.hpp"
#include "../numeric/units.hpp"

namespace flick { 
    begin_test_case(ocean_test_A) {
    using namespace material;
    ocean::configuration c;
    c.set<size_t>("n_angles",100);
    auto oce = ocean{c};
    size_t n1 = oce.material_ids().size();
    c.set<std::string>("mp_names",{"SD16_VF17","SD16_VF18"});
    c.set<double>("mp_concentrations",{1,1});
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
    c.set<double>("surface_depth_fraction",1);
    c.set<double>("cdom_440",10);
    auto oce = ocean{c};
    oce.set_wavelength(wl);
    oce.set_position({0,0,0});
    oce.set_direction({0,0,-1});
    double od1 = oce.absorption_optical_depth(distance);
    c.set<double>("surface_depth_fraction",0.5);
    oce = ocean{c};
    oce.set_wavelength(wl);
    oce.set_position({0,0,0});
    oce.set_direction({0,0,-1});
    double od2 = oce.absorption_optical_depth(distance);
    check_close(od1/od2,2,0.1_pct);
  } end_test_case()
}
