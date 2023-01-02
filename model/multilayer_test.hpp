#include "multilayer.hpp"
#include "../material/henyey_greenstein.hpp"
#include "../material/water/pure_water.hpp"

namespace flick {
  begin_test_case(multilayer_test) {
    using namespace flick;
  
    model::layer l = model::bottom_layer<coating::grey_lambert>(0,0);
    l.activate_receivers();

    model::plane_parallel_structure s(l);

    l = model::layer{thickness{1},"ocean"};
    absorption_coefficient a{0};
    scattering_coefficient b{0.0};
    asymmetry_factor g{0.0};
    
    l.fill<material::pure_water>();
    l.activate_receivers();
    s.add_on_top(l);
    l.name("atmosphere");
    l.fill<material::henyey_greenstein>(a,b,g);
    s.add_on_top(l);
    
    //std::cout << s; 

    emitter em{{0,0,0.0},100};
    em.set_direction<unidirectional>(unit_vector{0,0,1});
    s.transport_radiation(em,"ocean");
    receiver& r = s.outward_receiver("atmosphere");
    //std::cout << r.radiant_flux();
  } end_test_case()
}
