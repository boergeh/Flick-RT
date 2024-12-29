#include "marine_particles.hpp"

namespace flick {
  begin_test_case(marine_particles_test) {
    double massc = 1e-3;
    double wl = 550e-9;
    std::string material_name = "SD16_VF17";
    material::marine_particles mp(material_name);
    mp.set_wavelength(wl);
    mp.mass_concentration(massc);
    check_close(mp.scattering_coefficient(),0.8817,0.1);
    
    auto p_forward = mp.mueller_matrix(unit_vector{0.00,0}).value(0,0);
    auto unit_back = unit_vector{constants::pi,0};
    auto p_backward = mp.mueller_matrix(unit_back).value(0,0);
    
    check(p_forward/p_backward > 1e4);
    check_close(material::phase_function(mp).integral_4pi(),1,1_pct);
    check_close(material::phase_function(mp).asymmetry_factor(),0.943,5);

    double bleaching = 1;
    double scat_scale = 1;
    material::marine_particles mp_bleached(material_name,massc,scat_scale,bleaching);
    mp_bleached.set_wavelength(wl);
    check(mp_bleached.absorption_coefficient() < mp.absorption_coefficient());

    bleaching = 1e14;
    material::marine_particles mp_bleached2(material_name,massc,scat_scale,bleaching);
    check_small(mp_bleached2.absorption_coefficient());
  } end_test_case()
}
