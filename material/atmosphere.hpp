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
	
	add<double>("water_vapor", 17, R"(Vater vapor column thickness [m] at STP.)");
	
	add<double>("aerosol_od", 0, R"(Aerosol vertical total optical thickness at 550 nm)");
	
	add<double>("aerosol_ratio", 1, R"(Ratio of rural aerosol optical depth to total aerosol optical
depth. Set to '1' for rural aerosols only and '0' for urban aerosols only)");
	
	add<double>("relative_humidity", 0.5, R"(Relative humidity ratio of the atmosphere at the ground)");
	
	add<double>("cloud_liquid", 0, R"(Liquid equivalent cloud thickness [m]. Typical values are in the range
0 to 1e-4 m)");
	
	add<double>("snow_ice", 0, R"(Ice equivalent snow thickness [m], which is the thickness the snow
layer would have if it was compressed into a layer of ice with no air
between snow grains.)");
	
	add<double>("snow_radius", 100e-6, R"(Average snow grain radius [m].)");
	
	add<std::string>("gases", {"o3","o2","h2o"}, R"(Space-separated list of names of absorbing gases included in the
atmosphere, selected among 'o3', 'o2', 'h2o', 'no2', and 'co2')");
	
	add<std::string>("gas_spectral_region", "uv_vis", R"(Atmospheric spectral region for pre-calculated smoothed gas absorption
spectra. Select 'uv_vis' or 'uv_vis_toa', where 'uv_vis_toa' is
optimized for TOA radiation.)");
      }
    };
  private:
    basic_configuration c_;
    static const size_t n_cloud_base_ = 1;
  public:
    atmosphere(const basic_configuration& c=atmosphere::configuration())
      : mixture(angle_range(c.get<size_t>("n_angles")), height_grid(c)) {
      c_ = c;
      auto_update_iops(false);
      add_air();
      add_aerosols();
      add_clouds();
      add_snow();
      auto_update_iops(true);
    }
    static stdvector height_grid(const basic_configuration& c) {
      double epsilon = 1e-4; // m
      stdvector h = atmospheric_state(c.get<size_t>("n_heights")).height_grid();
      double cl = c.get<double>("cloud_liquid");
      double aod = c.get<double>("aerosol_od");
      if (cl > 0 or aod > 0) { 
	double h_cloud_top = h.at(n_cloud_base_+1);
	double h_cloud_base = h.at(n_cloud_base_);
	h.insert(h.begin()+n_cloud_base_,h_cloud_base-epsilon); 
	h.insert(h.begin()+n_cloud_base_+3,h_cloud_top+epsilon); 
      }
      double si = c.get<double>("snow_ice");     
      if (si > 0) { // layer and step layer for snow
	double snow_depth = 1;
	h.insert(h.begin()+1,snow_depth+epsilon); 
	h.insert(h.begin()+1,snow_depth); 
      }
      return h;
    }
  private:
    size_t cloud_base_above_step() {
      size_t n = n_cloud_base_ + 1;
      if (c_.get<double>("snow_ice") > 0) {
	return n + 2;
      }
      return n;
    }
    double layer_thickness(size_t n_base, size_t n_top) {
    	return heights().at(n_top)-heights().at(n_base);
    }
    void add_air() {
      double p = c_.get<double>("pressure");
      if (p > 0) {
	atmospheric_state s(c_.get<double>("temperature"),p);
	s.remove_all_gases();
	for (size_t i=0; i<c_.size<std::string>("gases"); i++) {
	  std::string gas = c_.get<std::string>("gases",i);
	  s.add_gas(gas);
	  if (gas=="o3")	  
	    s.scale_to_stp_thickness("o3",c_.get<double>("ozone"));
	  if (gas=="h2o")	  
	    s.scale_to_stp_thickness("h2o",c_.get<double>("water_vapor"));
	}	
	add_material<smooth_air>(s,c_.get<std::string>("gas_spectral_region"));
      }
      else
	add_material<vacuum>();
    }
    void add_aerosols() {
      double aod = c_.get<double>("aerosol_od");
      if (aod > 0) {
	size_t n_base = 0;
	size_t n_top = cloud_base_above_step()-1;     
	double ratio = c_.get<double>("aerosol_ratio");
	double rh = c_.get<double>("relative_humidity");
	double dh = layer_thickness(n_base,n_top);
	add_material<rural_aerosols>(dh, aod*ratio, rh);
	add_material<urban_aerosols>(dh, aod*(1-ratio), rh);
	set_range<rural_aerosols>(n_base,n_top);
	set_range<urban_aerosols>(n_base,n_top);
      }
    }
    void add_clouds() {
      double liquid_depth = c_.get<double>("cloud_liquid");
      if (liquid_depth > 0) {
	size_t n_base = cloud_base_above_step();
	size_t n_top = n_base + 1; 
	double radius = 15e-6;
	double dh = layer_thickness(n_base,n_top);
	double volume_fraction = liquid_depth/dh;  
	double mu = log(radius);
	double sigma = 0;
	using cloud = water_cloud<parameterized_monodispersed_mie>;
	add_material<cloud>(volume_fraction,mu,sigma);
	set_range<cloud>(n_base,n_top);
      }
    }
    void add_snow() {
      double ice_depth = c_.get<double>("snow_ice");
      if (ice_depth > 0) {
	size_t n_base = 0;
	size_t n_top = 1;
	double radius = c_.get<double>("snow_radius");
	double dh = heights().at(n_top)-heights().at(n_base);
	double volume_fraction = ice_depth/dh;  
	double mu = log(radius);
	double sigma = 0;
	using snow = ice_cloud<parameterized_monodispersed_mie>;
	add_material<snow>(volume_fraction,mu,sigma);
	set_range<snow>(n_base,n_top);
      }
    }
  };
}
}

#endif
