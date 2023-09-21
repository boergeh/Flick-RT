#include "mixture.hpp"
#include "../numeric/legendre/delta_fit.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;
  
  begin_test_case(mixture_test_A) {
    
  } end_test_case()
  
  begin_test_case(mixture_test) {
    double pi = constants::pi;
    double cloud_scat_coef = 1e-4;
    double zc_low = 100;
    double zc_high = 700;
    double toa = 120e3;
    stdvector heights = {0,100,700,800,2000,toa};
    stdvector angles = range(0,pi,100).linspace();

    auto sky =material::mixture(angles, heights);
    using simple_cloud = material::white_isotropic;
   
    sky.add_material<simple_cloud>(cloud_scat_coef);
    sky.set_range<simple_cloud>(zc_low,zc_high);
    sky.update_iops();
    
    double od_cl = sky.optical_depth(toa);
    double od_cl_bench = cloud_scat_coef*(zc_high-zc_low);
    check_close(od_cl,od_cl_bench);
    double p = sky.mueller_matrix(unit_vector{pi/2,0}).value(0,0);
    check_close(p,1/(4*pi),1e-7_pct);
    auto pf = material::phase_function(sky);
    check_close(delta_fit(pf,8).coefficients()[0]*4*pi,1,1e-7_pct);
    sky.add_material<material::rural_aerosols>();
    sky.update_iops();
    double od_ae = sky.get_material<material::rural_aerosols>().optical_depth(toa);
    check_close(sky.optical_depth(toa),od_ae+od_cl_bench);
    auto p1 = sky.mueller_matrix(unit_vector{0,0}).value(0,0);
    auto p2 = sky.mueller_matrix(unit_vector{pi,0}).value(0,0);
    check(p1 > 2/(4*pi));
    check(p2/p1 < 0.3);
    auto pf2 = material::phase_function(sky);
    //std::cout << "phase function terms: "<<std::setprecision(5) <<delta_fit(pf2,16).coefficients()*4*pi;
    check_close(delta_fit(pf2,16).coefficients()[0]*4*pi,1,0.05_pct);

  } end_test_case()
}
