#ifndef flick_material_atmosphere
#define flick_material_atmosphere

#include "gas/air.hpp"
#include "gas/atmospheric_state.hpp"
#include "aerosols/aerosols.hpp"
#include "spheres.hpp"
#include "mixture.hpp"

namespace flick {
namespace material {
  struct atmosphere : public mixture {
    struct config {
      size_t angles = 0;
      size_t heights = 8;
      double surface_temperature = 290;
      double surface_pressure = 1000e2;
      double ozone_stp_thickness = 0.003; 
      double aerosol_optical_depth_at_550_nm = 0.1;
      double aerosol_rural_to_urban_ratio = 1;
      double relative_humidity = 0.5;
      double cloud_liquid_thickness = 0.001;
      std::string wavelength_range = "uv_vis";
      std::vector<std::string> gases = {"o3","o2","h2o","no2"};
    };

    atmosphere(const atmosphere::config& c)
      : mixture(range(0,constants::pi,c.angles).linspace(),
		atmospheric_state(c.heights).height_grid()) {
      atmospheric_state s(c.surface_temperature,c.surface_pressure);
      s.remove_all_gases();
      for (size_t i=0; i<c.gases.size(); i++) {
	s.add_gas(c.gases[i]);
      }
      s.scale_to_stp_thickness("o3",c.ozone_stp_thickness);
      should_update_iops(false);
      set_material<smooth_air>(s,c.wavelength_range);
      double r = c.aerosol_rural_to_urban_ratio;
      set_material<rural_aerosols>(c.aerosol_optical_depth_at_550_nm*r,
       				   c.relative_humidity);
      set_material<urban_aerosols>(c.aerosol_optical_depth_at_550_nm*(1-r),
       				   c.relative_humidity);
      double base = 3e3;
      double top = 4e3;
      double radius = 5e-6;
      double volume_fraction = c.cloud_liquid_thickness/(top-base);  
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
