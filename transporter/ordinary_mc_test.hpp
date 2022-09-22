#include "ordinary_mc.hpp"
#include "../component/emitter.hpp"
#include "../material/henyey_greenstein.hpp"

namespace flick {
  begin_test_case(ordinary_mc_test_A) {
    semi_infinite_box box;
    box.name("box");
    box.move_by({0,0,1});
    box().outward_receiver().activate();
    emitter emitter{2};
    emitter.set_direction<unidirectional>(unit_vector{0,0,1});
    transporter::ordinary_mc omc(box);
    omc.transport_radiation(emitter,"box");
    check_close(omc.outward_receiver("box").radiant_flux(),2,1);
  } end_test_case()

  begin_test_case(ordinary_mc_test_B) {
    semi_infinite_box box;
    box.name("box");
    box.move_by({0,0,1});
    box().outward_receiver().activate();
    box().coat<coating::grey_lambert>(0,1);
    size_t np = 3;
    emitter emitter{np};
    emitter.set_direction<unidirectional>(unit_vector{0,0,1});
    transporter::ordinary_mc omc{box};
    omc.transport_radiation(emitter,"box");
    check_close(omc.outward_receiver("box").radiant_flux(),np,1);
  } end_test_case()
 
  begin_test_case(ordinary_mc_test_C) {
    semi_infinite_box layer;
    semi_infinite_box bottom;
    layer.name("layer");
    layer().outward_receiver().activate();
    bottom().inward_receiver().activate();
    layer.move_by({0,0,1});
    bottom().coat<coating::grey_lambert>(1,0);
    layer.insert(bottom);
    size_t np = 1;
    emitter emitter{{0,0,0.5},np};
    emitter.set_direction<unidirectional>(unit_vector{0,0,-1});
    transporter::ordinary_mc omc{layer};
    omc.transport_radiation(emitter,"layer");
    check_close(omc.outward_receiver("layer").radiant_flux(),np,1);
  } end_test_case()

  begin_test_case(ordinary_mc_test_D) {
    double r = 1;
    sphere s(r);
    s.name("s");
    s().outward_receiver().activate();

    size_t n = 2000;
    emitter emitter{n};
    emitter.set_direction<unidirectional>(unit_vector{0,0,1});
    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0};
    
    s().fill<material::henyey_greenstein>(a,b,g);
    transporter::ordinary_mc omc_a{s};
    omc_a.transport_radiation(emitter,"s");
    check_close(omc_a.outward_receiver("s").radiant_flux(),n,1,"a");
    check_close(omc_a.outward_receiver("s").mean_traveling_length(),r,1,"b");

    a = 1;
    s().fill<material::henyey_greenstein>(a,b,g);
    transporter::ordinary_mc omc_b{s};
    omc_b.transport_radiation(emitter,"s");
    check_close(omc_b.outward_receiver("s").radiant_flux(),n*exp(-r),0.1,"c");
    
    a = 0;
    b = 1;
    g = 0.5;
    s().fill<material::henyey_greenstein>(a,b,g);
    transporter::ordinary_mc omc_c{s};
    omc_c.transport_radiation(emitter,"s",g());
    check_close(omc_c.outward_receiver("s").radiant_flux(),n,9,"d");

    transporter::ordinary_mc omc_d{s};
    omc_d.transport_radiation(emitter,"s",g()-0.1);
    check_close(omc_d.outward_receiver("s").radiant_flux(),n,9,"e");
    
  } end_test_case()
  
  begin_test_case(ordinary_mc_test_E) {
    semi_infinite_box outer_box;
    semi_infinite_box inner_box;
    semi_infinite_box bottom_box;
    inner_box.name("inner_box");
    outer_box.name("outer_box");
    outer_box.move_by({0,0,3});
    inner_box.move_by({0,0,2});
    absorption_coefficient a{0};
    scattering_coefficient b{0};
    asymmetry_factor g{0};
    double real_n{9};
    inner_box().fill<material::henyey_greenstein>(a,b,g,real_n);
    inner_box().inward_receiver().activate();
    bottom_box().inward_receiver().activate();
    inner_box().outward_receiver().activate();
    outer_box.insert(inner_box);

    size_t n = 2000;
    emitter emitter{{0,0,2.5},n};
    emitter.set_direction<unidirectional>(unit_vector{0,0,-1});   
    transporter::ordinary_mc omc{outer_box};
    omc.transport_radiation(emitter,"outer_box");
    double reflectance = omc.outward_receiver("inner_box").radiant_flux() / n;
    double transmittance = omc.inward_receiver("inner_box").radiant_flux() / n;
    double r_benchmark = pow((1-real_n)/(1+real_n),2);
    double t_benchmark = 1-r_benchmark;
    check_close(reflectance,r_benchmark,5,"r");
    check_close(transmittance,t_benchmark,5,"t");

  } end_test_case()
}
