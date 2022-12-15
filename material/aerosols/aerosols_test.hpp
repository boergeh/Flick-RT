#include "aerosols.hpp"

namespace flick {
  begin_test_case(aerosols_test) {
    using namespace flick;
    material::rural_aerosols ra;
    ra.set(pose{{0,0,1},unit_vector{0,0,1}});
    ra.set(wavelength{500e-9});
    check_close(ra.scattering_coefficient(),1e-4,40);
    check_close(ra.absorption_coefficient(),1e-5,40);
    //std::cout << ra;
    
  } end_test_case()
}
