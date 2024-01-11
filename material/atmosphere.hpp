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
  class atmosphere : public mixture<pe_function> {
  public:
    struct configuration : public mixture::configuration {
      configuration() {
	add<double>("temperature", 290, R"(Atmosphere ground temperature [K])");
		    
	add<double>("pressure", 1000e2, R"(Atmosphere ground pressure [Pa])");
	
	add<double>("ozone", 0.003, R"(Ozone column thickness [m] at STP. Note that 100 DU = 0.001 m)");
	
	add<double>("aerosol_od", 0.01, R"(Aerosol vertical total optical thickness at 550 nm)");
	
	add<double>("aerosol_ratio", 1, R"(Ratio of rural aerosol optical depth to total aerosol optical
depth. Set to '1' for rural aerosols only and '0' for urban aerosols only)");
	
	add<double>("relative_humidity", 0.5, R"(Relative humidity ratio of the atmosphere at the ground)");
	
	add<double>("cloud_liquid", 0, R"(Liquid equivalent cloud thickness [m]. Typical values are in the range
0 to 1e-4 m)");
	
	add<std::string>("gases", {"o3","o2","h2o","no2"}, R"(Space-separated list of names of absorbing gases included in the
atmosphere, selected among 'o3', 'o2', 'h2o', and 'no2')");
      }
    };
  private:
    basic_configuration c_;
  public:
    atmosphere(const basic_configuration& c=atmosphere::configuration())
      : mixture(angle_range(c.get<size_t>("n_angles")),
		atmospheric_state(c.get<size_t>("n_heights")).height_grid()) {
      c_ = c;
      add_air();
      add_aerosols();
      add_clouds();
      //update_iops();
    }
  private:
    void add_air() {
      atmospheric_state s(c_.get<double>("temperature"),c_.get<double>("pressure"));
      s.remove_all_gases();
      for (size_t i=0; i<c_.size<std::string>("gases"); i++) {
	s.add_gas(c_.get<std::string>("gases",i));
      }
      s.scale_to_stp_thickness("o3",c_.get<double>("ozone"));
      add_material<smooth_air>(s,"uv_vis");
    }
    void add_aerosols() {
      double ratio = c_.get<double>("aerosol_ratio");
      double aod = c_.get<double>("aerosol_od");
      double rh = c_.get<double>("relative_humidity");
      add_material<rural_aerosols>(aod*ratio, rh);
      add_material<urban_aerosols>(aod*(1-ratio), rh);
    }
    void add_clouds() {
      size_t n_base = 1;
      size_t n_top = 2;
      double radius = 15e-6;
      double dh = heights().at(n_top)-heights().at(n_base);
      double volume_fraction = c_.get<double>("cloud_liquid")/dh;  
      double mu = log(radius);
      double sigma = 0;
      using cloud = water_cloud<parameterized_monodispersed_mie>;
      set_range<cloud>(n_base,n_top);
      add_material<cloud>(volume_fraction,mu,sigma);
    }
  };
}
}

#endif
