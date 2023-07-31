#include "atmosphere.hpp"
#include "../numeric/units.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    using namespace units;
    using namespace material;
    atmospheric_state state(290_K, 1000_hPa, 8);
    state.remove_gas("co2");
    state.remove_gas("h2o");
    state.remove_gas("o2");

    atmosphere atm{2, state.height_grid()};
    atm.set_material<air>(state);
    atm.set_wavelength(350_nm);
    
    //atm.get_material<air>().set_wavelength(300_nm);
    //auto& m = atm.get_material<rural_aerosols>();
    //m.set_wavelength(350_nm);
    atm.update_iops();
    std::cout << atm.get_material<rural_aerosols>().scattering_optical_depth(100e3);
    
  } end_test_case()
}
