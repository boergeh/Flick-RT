#include "mixture.hpp"
#include "atmosphere.hpp"
#include "ocean.hpp"
#include "../numeric/legendre/delta_fit.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;
  
  begin_test_case(mixture_test_A) {
    double pi = constants::pi;
    double cloud_scat_coef = 1e-4;
    size_t n_low = 1;
    size_t n_high = 2;
    double toa = 120e3;
    stdvector heights = {0,100,700,800,2000,toa};
    stdvector angles = range(0,pi,100).linspace();

    auto sky =material::mixture(angles, heights);
    using simple_cloud = material::white_isotropic;
   
    sky.add_material<simple_cloud>(cloud_scat_coef);
    sky.set_range<simple_cloud>(n_low, n_high);
    double od_cl = sky.optical_depth(toa);
    double od_cl_bench = cloud_scat_coef*(heights[n_high]-heights[n_low]);
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
    check_close(delta_fit(pf2,16).coefficients()[0]*4*pi,1,0.07_pct);
  } end_test_case()
  
    begin_test_case(mixture_test_B) {
    using namespace material;
    material::mixture m({0,1,3.14},{-100,-1e-6,0,10e3,50e3,100e3});
    //m.should_update_iops(false);
    m.add_material<ocean>();
    m.add_material<atmosphere>();
    m.set_range<ocean>(0,1);
    m.set_range<atmosphere>(2,5);
    //m.should_update_iops(true);
    //m.update_iops();
  } end_test_case()
}
