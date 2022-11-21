#include "profile.hpp"

namespace flick {
  begin_test_case(profile_test) {
    using namespace constants;
    using namespace flick;
    profile::constant c{1,pose{{0,0,0},{pi/4,0}}};
    check_close(c.optical_depth(2),2,1e-12);
    check_close(c.distance(3),3,1e-12);

    double K = 1/9.0;

    // surface scattering coefficient 1 km^-1;
    // point two at 10 km
    pe_function f{{0,10},{1, exp(-K*10)}};
    profile::z_varying p_vertical{f,pose{{0,0,0},{0,0,1}}};
    double z_max = 20;
    double tau_zenith = (1-exp(-K*z_max))/K;
    check_close(p_vertical.optical_depth(z_max),tau_zenith,1e-9,"a");

    // See Wikipedia air mass
    double ang_slant = 88.57/180*pi;
    double tau_slant = tau_zenith * 40;
    profile::z_varying p_slant_up{f,pose{{0,0,0},{ang_slant,0}}};
    double d = p_slant_up.distance(tau_slant);
    check_close(p_slant_up.optical_depth(d),tau_slant,1e-9,"b");

    double z_toa = d*sin(pi/2-ang_slant);
    profile::z_varying p_slant_down{f,pose{{0,0,z_toa},{pi-ang_slant,0}}};
    check_close(p_slant_down.distance(tau_slant),d,1e-9,"c");
    check_close(p_slant_down.optical_depth(d),tau_slant,1e-9,"d");

    profile::z_varying p_horizontal{f,pose{{0,0,0},{pi/2,0}}};
    check_close(p_horizontal.optical_depth(1),1,1e-9,"e");
    check_close(p_horizontal.distance(1),1,1e-9,"e");
    
    double tau_horizon = tau_zenith * 35;
    //    profile::z_varying zv{f,pose{{0,0,0},unit_vector{pi/4,0}}};
  } end_test_case()
}
