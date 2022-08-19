#include "mc_basic.hpp"

namespace flick {  

  /*
  geometry::volume<content> world() {
    using semi_infinite_box = geometry::semi_infinite_box<content>;
    semi_infinite_box space;
    semi_infinite_box atmosphere;
    semi_infinite_box bottom;
  
    space.name("space");
    atmosphere.name("atmosphere");
    bottom.name("bottom");
    
    space.content().inward_receiver().activate();
    space.content().outward_receiver().activate();
    atmosphere.content().inward_receiver().activate();
    atmosphere.content().outward_receiver().activate();
    bottom.content().inward_receiver().activate();
    bottom.content().outward_receiver().activate();
    bottom.content().set_coating<coating::grey_lambert>(1,0);
  
    space.move_by({0,0,10});
    atmosphere.move_by({0,0,1});
    atmosphere.insert(bottom);
    space.insert(atmosphere);
  
    return space;
  }
  */
  
  begin_test_case(mc_basic_test) {
    //geometry::volume<content> w = world();
    //mc_basic mcb{w};
    //mcb.run();
    /*
    geometry::plane_parallel::henyey_greenstein_slab hgs;
    hgs.n_radiation_packages(5);   
    check_small(hgs.hemispherical_reflectance(),1e-12);
    check_close(hgs.hemispherical_transmittance(),1,1e-12);

    hgs.bottom_albedo(1);
    check_close(hgs.hemispherical_reflectance(),1,1e-12);
    check_small(hgs.hemispherical_transmittance(),1e-12);
    */
    } end_test_case()
}
