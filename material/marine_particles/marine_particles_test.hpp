#include "marine_particles.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(marine_particles_test) {
    using namespace units;
    material::marine_particles mp("SD16_VF17");
    mp.set_wavelength(550e-9);
    mp.mass_concentration(1e-3);
    check_close(mp.scattering_coefficient(),0.8817,0.1);
    auto p_forward = mp.mueller_matrix(unit_vector{0.00,0}).value(0,0);
    auto unit_back = unit_vector{constants::pi,0};
    auto p_backward = mp.mueller_matrix(unit_back).value(0,0);
    check(p_forward/p_backward > 1e4);
    check_close(material::phase_function(mp).integral_4pi(),1,1.6);
    check_close(material::phase_function(mp).asymmetry_factor(),0.943,0.1);  
  } end_test_case()
}
