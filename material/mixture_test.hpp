#include "mixture.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;
  
  begin_test_case(mixture_test) {
    double cloud_scat_coef = 1e-4;
    double zc_low = 100;
    double zc_high = 700;
    double toa = 120e3;
    stdvector heights = {0,100,700,800,2000,toa};
    stdvector angles = {0,3.14};

    using namespace material;
    mixture sky{angles, heights};
    using simple_cloud = white_isotropic;

    sky.add_material<simple_cloud>(cloud_scat_coef);
    sky.set_scaling_factor<simple_cloud>(layer_only(zc_low,zc_high));
    sky.update_iops();
    
    double od_cl = sky.optical_depth(toa);
    double od_cl_bench = cloud_scat_coef*(zc_high-zc_low);
    check_close(od_cl,od_cl_bench);

    sky.add_material<rural_aerosols>();
    sky.update_iops();
    double od_ae = sky.get_material<rural_aerosols>().optical_depth(toa);
    check_close(sky.optical_depth(toa),od_ae+od_cl_bench);  
  } end_test_case()
}
