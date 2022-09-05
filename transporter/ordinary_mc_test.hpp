#include "ordinary_mc.hpp"
#include "../component/emitter.hpp"
#include "../material/henyey_greenstein.hpp"

namespace flick {
  begin_test_case(ordinary_mc_test_A) {
    semi_infinite_box box;
    box.move_by({0,0,1});
    box().outward_receiver().activate();
    box().coat<coating::grey_lambert>(0,1);
    size_t np = 1;
    emitter emitter{np};
    emitter.set_direction<unidirectional>(unit_vector{0,0,1});
    transporter::ordinary_mc(box,box,emitter).run();
    check_close(box().outward_receiver().radiant_flux(),np,1);
  } end_test_case()
 
  begin_test_case(ordinary_mc_test_B) {
    semi_infinite_box layer;
    semi_infinite_box bottom;
    layer().outward_receiver().activate();
    bottom().inward_receiver().activate();
    layer.move_by({0,0,1});
    bottom().coat<coating::grey_lambert>(1,0);
    layer.insert(bottom);
    size_t np = 1;
    emitter emitter{{0,0,0.5},np};
    emitter.set_direction<unidirectional>(unit_vector{0,0,-1});
    transporter::ordinary_mc(layer,layer,emitter).run();
    check_close(layer().outward_receiver().radiant_flux(),np,1);
  } end_test_case()
  
  begin_test_case(ordinary_mc_test_C) {
    double r = 1;
    sphere s(r);
    s().outward_receiver().activate();

    size_t n = 10;
    emitter emitter{n};

    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0};
    s().fill<material::henyey_greenstein>(a,b,g);
    transporter::ordinary_mc(s,s,emitter).run();
    check_close(s().outward_receiver().radiant_flux(),n,1);
    check_close(s().outward_receiver().mean_traveling_length(),r,1,"a");

    a = 1;
    emitter.add_packages(n);
    s().outward_receiver().clear();
    s().fill<material::henyey_greenstein>(a,b,g);
    transporter::ordinary_mc(s,s,emitter).run();
    check_close(s().outward_receiver().radiant_flux(),n*exp(-r),1,"b");

    a = 0;
    b = 1;
    g = 0.5;
    emitter.add_packages(n);
    s().outward_receiver().clear();
    s().fill<material::henyey_greenstein>(a,b,g);
    transporter::ordinary_mc(s,s,emitter).run(1,0.5);
    check_close(s().outward_receiver().radiant_flux(),n,1,"c");

    n *= 100; 
    emitter.add_packages(n);
    s().outward_receiver().clear();
    s().fill<material::henyey_greenstein>(a,b,g);
    transporter::ordinary_mc(s,s,emitter).run(1,0.6);
    check_close(s().outward_receiver().radiant_flux(),n,5,"diff.sampl.g");
      
  } end_test_case()
}
