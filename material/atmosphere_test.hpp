#include "atmosphere.hpp"
#include "../numeric/units.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    using namespace units;
    using namespace material;
    //atmospheric_state state(290_K, 1000_hPa, 8);
    //state.remove_gas("co2");
    //state.remove_gas("h2o");
    // state.remove_gas("o2");
    // state.remove_gas("o3");
    // std::cout << "ppm: " << std::setprecision(4)<<state.fraction("co2")*1e6 << std::endl;
    //state.scale_to_fraction("co2",180e-6);
    //std::cout << "ppm: " << std::setprecision(4)<<state.fraction("co2")*1e6 << std::endl;
    
    //atmosphere atm{0, state.height_grid()};
    //atm.set_material<air>(state);

  } end_test_case()
}
