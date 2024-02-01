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
}
