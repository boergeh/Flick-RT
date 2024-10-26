#include "air.hpp"
#include "../../numeric/units.hpp"

namespace flick {
  begin_test_case(air_test_o3)
  // Van Weele M, et al., Journal of Geophysical Research:
  // Atmospheres, 105(D4), pp.4915-4925, Table 5.
  {
    using namespace units;
    atmospheric_state state(290_K, 1019_hPa, 8);
    state.remove_all_gases();
    state.add_gas("o3");
    state.scale_to_stp_thickness("o3", 2.95_mm);

    material::hitran_air air(state);
    air.set_wavelength(310_nm);
    air.set_position({0,0,17_m});
    air.set_direction({0,0,1});
    check_close(air.scattering_optical_depth(120_km), 1.059, 0.2_pct);
    check_close(air.absorption_optical_depth(120_km), 0.687, 0.2_pct);

    material::smooth_air air_uv(state,"uv");
    air_uv.set_wavelength(310_nm);
    air_uv.set_position({0,0,17_m});
    air_uv.set_direction({0,0,1});
    check_close(air_uv.scattering_optical_depth(120_km), 1.059, 1_pct);
    check_close(air_uv.absorption_optical_depth(120_km), 0.687, 1_pct);

    air_uv.set_wavelength(3000_nm);
    check(air_uv.scattering_optical_depth(120_km)<0.001);
    
  } end_test_case()

  begin_test_case(air_test_o2)
  // From Qilong Min, replotted in C.F. Bohern and E.E. Clothiaux,
  // Figure 6.14, p 321.
  {
    using namespace units;
    bool plot_spectrum = false;
    atmospheric_state state(300_K, 1000_hPa, 8);
    state.remove_all_gases();
    state.add_gas("o2");

    material::hitran_air air(state);
    air.set_wavelength(1./13122*1e-2);
    check_close(air.absorption_optical_depth(100_km), 0.037, 10_pct);
    air.set_wavelength(1./12965*1e-2);
    check_close(air.absorption_optical_depth(100_km), 0.11, 85_pct);

    if (plot_spectrum) {
      std::ofstream of("./tmp.txt");
      double wl2 = 1./12870*1e-2;
      double wl1 = 1./13170*1e-2;
      for (auto& wl:range(wl1,wl2,50000).logspace()) {
	air.set_wavelength(wl);
	of << std::setprecision(7);
	of << wl << " " << air.absorption_optical_depth(100_km) << std::endl;
      }
      of.close();
    }
  } end_test_case()
  

}
