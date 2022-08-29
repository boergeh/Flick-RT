#include "single_layer_slab.hpp"

namespace flick {  
  begin_test_case(single_layer_slab_test) {
    using namespace flick;
    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0.5};
    thickness h{1};
    //std::cout << g;
    model::single_layer_slab slab{h};
    slab.fill<material::henyey_greenstein>(a,b,g);    
    slab.set(bottom_albedo{0});
    slab.set(incidence_angle{0});
    //slab.set(precision(1));
    slab.set(number_of_packages{100});

    check_small(slab.hemispherical_reflectance(),1e-12,"a");

    slab.set(bottom_albedo{1});
    check_close(slab.hemispherical_reflectance(),1,1e-12,"b");

    slab.set(bottom_albedo{0.5});
    check_close(slab.hemispherical_reflectance(),0.5,50,"c");

    a = 1;
    slab.fill<material::henyey_greenstein>(a,b,g);
    slab.set(bottom_albedo(0));
    check_close(slab.hemispherical_transmittance(),exp(-h()),1,"d");

    a = 0;
    b = 1;
    g = 0.5;
    slab.fill<material::henyey_greenstein>(a,b,g);
    slab.set(bottom_albedo{0});
    slab.set(incidence_angle{constants::pi/2*0.4});
    slab.hemispherical_reflectance();
    check_fast(0.1,"e");
    /*
    std::cout << slab.directional_reflectance(polar_angle{constants::pi/2},
					      azimuth_angle{0},
					      unit_interval{0.5});
    */
  } end_test_case()
}
