#include "iop_profile.hpp"

namespace flick {
  begin_test_case(iop_profile_test) {
    using namespace constants;
    using namespace flick;
    constant_iop c{1};
    check_close(c.optical_depth(2),2,1e-12);
    check_close(c.distance(3),3,1e-12);
    
    double K = 1/9.0;

    // surface scattering coefficient 1 km^-1;
    // point two at 10 km
    pe_function f = {{0,10},{1, exp(-K*10)}};
    pose start{{0,0,0},{0,0,1}};
    iop_z_profile atmlike{f};
    double z_max = 20;
    double tau_zenith = (1-exp(-K*z_max))/K;
    check_close(atmlike.optical_depth(start,z_max),tau_zenith,1e-9,"a");

    // See Wikipedia air mass
    double ang_slant = 88.57/180*pi;
    double tau_slant = tau_zenith * 40;
    start = pose{{0,0,0},{ang_slant,0}};
    double d = atmlike.distance(start,tau_slant);
    check_close(atmlike.optical_depth(start,d),tau_slant,1e-9,"b");

    double z_toa = d*sin(pi/2-ang_slant);
    check_close(atmlike.distance(start,tau_slant),d,1e-9,"c");
    check_close(atmlike.optical_depth(start,d),tau_slant,1e-9,"d");

    start = pose{{0,0,0},{pi/2,0}};
    check_close(atmlike.optical_depth(start,1),1,1e-9,"e");
    check_close(atmlike.distance(start,1),1,1e-9,"e");
    
    //double tau_horizon = tau_zenith * 35;
    //    profile::z_varying zv{f,pose{{0,0,0},unit_vector{pi/4,0}}};
  } end_test_case()
}
