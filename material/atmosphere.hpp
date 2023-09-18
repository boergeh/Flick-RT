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
    struct config : basic_configuration {
      config() {
	add<size_t>("angles", 0,
		    "Number of computed phase function angles");
	add<size_t>("heights", 8,
		    "Number of computed vertical gas profile height points");
	add<double>("temperature", 290,
		    "[K] Bottom of atmosphere ground temperature");
	add<double>("pressure", 1000e2,
		    "[Pa] Bottom of atmosphere ground pressure");
	add<double>("ozone", 0.003,
		    "[m] Ozone column thickness at stp, 100 DU = 0.001 m");
	add<double>("aerosol_od", 0.1,
		    "Vertical aerosol optical depth at 550 nm");
	add<double>("aerosol_ratio", 1,
		    "Rural to urban aerosol concentration ratio");
	add<double>("relative_humidity", 1,
		    "Surface relative humidity (for aerosols)");
	add<double>("cloud_liquid", 0.0001,
		    "[m] Cloud liquid thickness");
	add<std::string>("wavelength_range", "uv_vis",
	  "Name of smoothed gas absorption spectra with a give wavelength range");
	add<std::string>("gases", {"o3","o2","h2o","no2"},
			 "Name of gases to include");
      }
    };
   
    atmosphere(const atmosphere::config& c)
      : mixture(range(0,constants::pi,c.get<size_t>("angles")).linspace(),
		atmospheric_state(c.get<size_t>("heights")).height_grid()) {
      atmospheric_state s(c.get<double>("temperature"),c.get<double>("pressure"));
      s.remove_all_gases();
      for (size_t i=0; i<c.size<std::string>("gases"); i++) {
	s.add_gas(c.get<std::string>("gases",i));
      }
      s.scale_to_stp_thickness("o3",c.get<double>("ozone"));
      should_update_iops(false);
      add_material<smooth_air>(s,c.get<std::string>("wavelength_range"));
      double ratio = c.get<double>("aerosol_ratio");
      double aod = c.get<double>("aerosol_od");
      double rh = c.get<double>("relative_humidity");
      add_material<rural_aerosols>(aod*ratio, rh);
      add_material<urban_aerosols>(aod*(1-ratio), rh);
      double base = 3e3;
      double top = 4e3;
      double radius = 5e-6;
      double volume_fraction = c.get<double>("cloud_liquid")/(top-base);  
      double mu = log(radius);
      double sigma = 0;
      using cloud = water_cloud<parameterized_monodispersed_mie>;
      pe_function profile = pe_function{{base,top},{1,1}};
      set_scaling_factor<cloud>(profile);
      add_material<cloud>(volume_fraction,mu,sigma);
      should_update_iops(true);
      update_iops();
    }    
  };
}
}

#endif
