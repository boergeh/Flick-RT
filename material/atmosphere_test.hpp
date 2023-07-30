#include "atmosphere.hpp"
#include "../numeric/units.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    using namespace units;
    atmospheric_state state(290_K, 1000_hPa, 14);
    state.remove_gas("co2");
    state.remove_gas("h2o");
    state.remove_gas("o2");
    material::atmosphere atm{2, state.height_grid()};
    atm.set_material<material::air>(state);
    atm.get_material<material::air>().set_wavelength(300_nm);
    auto& m =  atm.get_material<material::rural_aerosols>();
    m.set_wavelength(350_nm);
    atm.update_iops();
    //std::cout << atm.get_material<material::rural_aerosols>().wavelength();
    
  } end_test_case()
}
