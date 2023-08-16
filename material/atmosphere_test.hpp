#include "atmosphere.hpp"
#include "../numeric/units.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    using namespace units;
    using namespace material;
    atmospheric_state state(290_K, 1000_hPa, 8);
    state.remove_gas("co2");
    //state.remove_gas("h2o");
    //state.remove_gas("o2");
    //state.remove_gas("o3");

    atmosphere atm{0, state.height_grid()};
    atm.set_material<air>(state);

    std::ofstream of("./tmp.txt");
    for (auto& wl:range(0.280e-6,2.5e-6,500).logspace()) {
      atm.set_wavelength(wl);
    
      //atm.get_material<air>().set_wavelength(300_nm);
      //auto& m = atm.get_material<rural_aerosols>();
      //m.set_wavelength(350_nm);
      atm.update_iops();
      of << std::setprecision(7);
      of << wl << " " << atm.absorption_optical_depth(100_km) << std::endl;
      //of << wl << " " << atm.scattering_optical_depth(100_km) << std::endl;
      //std::cout << atm.get_material<rural_aerosols>().scattering_optical_depth(100e3);
    }
    of.close();
  } end_test_case()
}
