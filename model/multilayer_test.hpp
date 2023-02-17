#include "multilayer.hpp"
#include "../material/henyey_greenstein.hpp"

namespace flick {
  begin_test_case(multilayer_test) {
    using namespace flick;
  
    model::layer l = model::bottom_layer<coating::grey_lambert>(0,0);
    l.activate_receivers();

    model::plane_parallel_structure s(l);

    l = model::layer{thickness{1},"hg1"};
    absorption_coefficient a{0};
    scattering_coefficient b{1};
    asymmetry_factor g{0};
    
    l.fill<material::henyey_greenstein>(a,b,g);
    l.activate_receivers();
    s.add_on_top(l);

    l.name("hg2");
    l.fill<material::henyey_greenstein>(a,b,g);
    s.add_on_top(l);
    
    //std::cout << s; 

    emitter em{{0,0,0},100};
    em.set_direction<unidirectional>(unit_vector{0,0,1});
    s.transport_radiation(em,"hg1");
    check(s.outward_receiver("hg2").radiant_flux()>0);
    //receiver& r = s.outward_receiver("hg2");
    //std::cout <<"\n" <<r.polar_angle_distribution(10,10);
  } end_test_case()
}
