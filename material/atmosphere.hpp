#ifndef flick_material_atmosphere
#define flick_material_atmosphere

#include "gas/air.hpp"
#include "gas/atmospheric_state.hpp"
#include "aerosols/aerosols.hpp"
#include "spheres.hpp"
#include "mixture.hpp"
#include "../environment/configuration.hpp"

namespace flick {
namespace material {
  struct atmosphere : public mixture {
    enum {angles, heights, temperature, pressure, ozone, aerosol_od, aerosol_ratio,
      relative_humidity, cloud_liquid, wavelength_range, gases};
    struct config : configuration {
      config() {
	set<size_t>(angles, 0, "Number of computed phase function angles");
	set<size_t>(heights, 8, "Number of computed vertical gas profile height points");
	set<double>(temperature, 290, "[K] Bottom of atmosphere ground temperature");
	set<double>(pressure, 1000e2, "[Pa] Bottom of atmosphere ground pressure");
	set<double>(ozone, 0.003, "[m] Ozone column thickness at stp, 100 DU = 0.001 m");
	set<double>(aerosol_od, 0.1, "Vertical aerosol optical depth at 550 nm");
	set<double>(aerosol_ratio, 1, "Rural to urban aerosol concentration ratio");
	set<double>(relative_humidity, 1, "Surface relative humidity (for aerosols)");
	set<double>(cloud_liquid, 0.0001, "[m] Cloud liquid thickness");
	set<std::string>(wavelength_range, "uv_vis", "Name of smoothed gas absorption spectra with a give wavelength range");
	set<std::string>(gases, {"o3","o2","h2o","no2"}, "Name of gases to include");
      }
    };
    /*
    std::map<std::string,int> parameter_names;
    struct config2 : configuration {
      config2() {
	add<size_t>("angles", 0, "Number of computed phase function angles");
	add_configuration(ocean::config()); //put in add_material..
      }
    };
    // should be add, not set?
    */
    atmosphere(const atmosphere::config& c)
      : mixture(range(0,constants::pi,c.get<size_t>(angles)).linspace(),
		atmospheric_state(c.get<size_t>(heights)).height_grid()) {
      atmospheric_state s(c.get<double>(temperature),c.get<double>(pressure));
      s.remove_all_gases();
      for (size_t i=0; i<c.size<std::string>(gases); i++) {
	s.add_gas(c.get<std::string>(gases,i));
      }
      s.scale_to_stp_thickness("o3",c.get<double>(ozone));
      should_update_iops(false);
      set_material<smooth_air>(s,c.get<std::string>(wavelength_range));
      double ratio = c.get<double>(aerosol_ratio);
      double aod = c.get<double>(aerosol_od);
      double rh = c.get<double>(relative_humidity);
      set_material<rural_aerosols>(aod*ratio, rh);
      set_material<urban_aerosols>(aod*(1-ratio), rh);
      double base = 3e3;
      double top = 4e3;
      double radius = 5e-6;
      double volume_fraction = c.get<double>(cloud_liquid)/(top-base);  
      double mu = log(radius);
      double sigma = 0;
      using cloud = water_cloud<parameterized_monodispersed_mie>;
      pe_function profile = pe_function{{base,top},{1,1}};
      set_scaling_factor<cloud>(profile);
      set_material<cloud>(volume_fraction,mu,sigma);
      should_update_iops(true);
      update_iops();
    }    
  };
}
}

#endif
