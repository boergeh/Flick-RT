#include "absorption_smoother.hpp"
#include "../../numeric/units.hpp"
#include "../../numeric/distribution.hpp"

namespace flick {
  begin_test_case(absorption_smoother_test)
  {
    using namespace units;
    struct source {
      double value(double wl) {
	return 1;
      }
    };
    atmospheric_state state(290_K, 1019_hPa, 8);

    distribution::normal d(300_nm, 5_nm);

    auto as = absorption_smoother(source(), state, "o3");

    stdvector center_wls = {280e-9};
    //std::cout << gaussian_smooth_cross_section(as, center_wls, 0.01, 0.1) << std::endl;
    
  } end_test_case()
}
