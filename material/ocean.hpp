#ifndef flick_material_ocean
#define flick_material_ocean

#include "water/pure_water.hpp"
#include "water/cdom.hpp"
#include "water/phytoplankton.hpp"
#include "water/nap.hpp"
#include "marine_particles/marine_particles.hpp"
#include "marine_cdom/marine_cdom.hpp"
#include "mixture.hpp"
#include "spheres.hpp"
#include "../environment/configuration.hpp"

namespace flick {
namespace material {
  struct ocean : public mixture<pl_function> {
    struct configuration : public mixture::configuration {
      configuration() {
	add<double>("bottom_depth", 200, R"(Total depth of the water column [m])");
	
	add<double>("concentration_relative_depths", {0,1}, R"(Space-separated list of depth fractions, df, from surface (df=0) to bottom (df=1)
defining a material scaling factor profile.)");
	
	add<double>("concentration_scaling_factors", {1,1}, R"(Space-separated list of scaling factors which scale the concentration
of all ocean materials except pure water.)");
	
	add<double>("cdom_440", 0.0, R"(CDOM absorption coefficient at 440 nm [1/m])");

	add<double>("cdom_slope", 0.017, R"(Slope of the CDOM absorption spectrum [1/nm]. Note the exception from
the SI unit convention)");
		    	
	add<double>("chl_concentration", 0, R"(Chlorophyll concentration in the water column [kg/m^3]. A
concentration of e.g., 10.0 mg/m^3 may be written as 10.0e-6 kg/m^3
for clarity)");

	add<double>("nap_concentration", 0, R"(Dry mass concentration of nonalgal particles in the water column
[kg/m^3]. A concentration of e.g., 10.0 g/m^3 may be written as
10.0e-3 kg/m^3 for clarity)");
	
	add<double>("bubble_volume_fraction", 0, R"(Bubble volume fraction in the water column [unitless])");

	add<double>("water_temperature", 290, R"(Temperature in the water column [K])");

	add<double>("water_salinity", 30, R"(Salinity of the water column [PSU])");

	add<std::string>("mp_names", "SD16_VF17", R"(Space-separated list of names of measured marine particles
with inherent optical properties tabulated in separate ASCII files stored
in the Flick directory material/marine_particles/iop_tables)");
	
	add<double>("mp_concentrations", 0, R"(Space-separated list of dry mass concentrations [kg/m^3] of measured
marine particles with inherent optical properties tabulated
in separated ASCII files in the Flick directory
material/marine_particles/iop_table, one concentration value for each
material given in mp_names. Note that a concentration of e.g., 10.0
g/m^3 may be written as 10.0e-3 kg/m^3 for clarity)");
	
	add<double>("mp_scattering_scaling_factors", 1, R"(Space-separated list of scaling factors [unitless] for manual scaling of the
scattering coefficient marine particles, one scaling factor for each listed marine particles name.)");
	
	add<double>("mp_bleaching_factors", 0, R"(Space-separated list of factors [unitless] for degree of particle
bleaching, one scaling factor for each listed marine particles name. 0
gives full absorption and 1 gives absorption after adding a bleaching
chemical. A factor larger than one will reduce the absorption beyond
the bleached values.)");
	
	add<std::string>("mcdom_names", "HF22_D001", R"(Space-separated list of names of measured marine CDOM with absorption
coefficients tabulated in separate ASCII files stored in the Flick
directory material/marine_cdom/iop_tables)");
	
	add<double>("mcdom_scaling_factors", 0, R"(Space-separated list of scaling factor for measured marine CDOM
absorption coefficients listed in separated ASCII files in the Flick
directory material/marine_cdom/iop_table, one concentration value for
each CDOM spectra given in mcdom_names.)");
      }
    };
  private:
    basic_configuration c_;
  public:
    ocean(const basic_configuration& c=ocean::configuration())
      : mixture(angle_range(c.get<size_t>("n_angles")),height_grid(c)) {
      c_ = c;    
      auto_update_iops(false);
      add_pure_water();
      add_cdom();
      add_phytoplankton();
      add_nap();
      add_bubbles();
      add_marine_particles();
      add_marine_cdom(); 
      auto_update_iops(true);
    }
    static stdvector height_grid(const basic_configuration& c) {
      double epsilon = 1e-6;      
      double depth = c.get<double>("bottom_depth");
      stdvector z = absolute_depth(c);
      stdvector h = {-depth};
      for (size_t i = 0; i < z.size(); ++i) {
	if (z[i] < -2*epsilon and z[i] > -depth+epsilon)
	  h.push_back(z[i]);
      }
      h.push_back(-epsilon);
      return h;
    }
  private:
    void add_pure_water() {
      add_material<pure_water>(c_.get<double>("water_salinity"),
			       c_.get<double>("water_temperature"));
    }
    void add_cdom() {
      double a440 = c_.get<double>("cdom_440");
      if (a440 > 0) {
	auto m = std::make_shared<cdom>(a440, c_.get<double>("cdom_slope"));
	add_profile(m, "cdom");
      }
    }
    void add_profile(const std::shared_ptr<base>& m, const std::string& name) {
      stdvector factor = c_.get_vector<double>("concentration_scaling_factors");
      std::reverse(factor.begin(),factor.end());
      stdvector z = absolute_depth(c_);
      add_material(make_scaled_z_profile<pl_function>(m,z,factor),name);
    }
    static stdvector absolute_depth(const basic_configuration& c) {
      double bottom_depth = c.get<double>("bottom_depth");
      stdvector relative_depth = c.get_vector<double>("concentration_relative_depths");
      std::reverse(relative_depth.begin(),relative_depth.end());
      return (-1)*relative_depth*bottom_depth;
    }
    void add_phytoplankton() {
      double chl = c_.get<double>("chl_concentration");
      if (chl > 0) {
	add_profile(std::make_shared<phytoplankton>(chl),"phytoplankton");
      }
    }
    void add_nap() {
      double con = c_.get<double>("nap_concentration");
      if (con > 0) {
	add_profile(std::make_shared<nap>(con),"nap");
      }
    }
    void add_bubbles() {
      double volume_fraction = c_.get<double>("bubble_volume_fraction");
      if (volume_fraction > 0) {
	double effective_radius = 50e-6;
	double mu = log(effective_radius);
	double sigma = 0;
	using bubbles = bubbles_in_water<parameterized_monodispersed_mie>;
	auto b = bubbles(volume_fraction,mu,sigma);
	add_profile(std::make_shared<bubbles>(volume_fraction,mu,sigma),"bubbles");
      }
    }
    void add_marine_particles() {
      std::vector<std::string> names = c_.get_vector<std::string>("mp_names");
      std::vector<double> concentrations = c_.get_vector<double>("mp_concentrations");
      std::vector<double> scattering_scaling_factors = c_.get_vector<double>("mp_scattering_scaling_factors");
      std::vector<double> bleaching_factors = c_.get_vector<double>("mp_bleaching_factors");
      //ensure(names.size() == concentrations.size() and names.size() ==
      //     scattering_scaling_factors.size() and names.size() ==
      //     bleaching_factors.size(), "marine particles");
      for (size_t i = 0; i < names.size(); i++) {
	if (concentrations.at(i) > 0) {
	  auto m = std::make_shared<marine_particles>(names.at(i),
						      at_or_last(concentrations,i),
						      at_or_last(scattering_scaling_factors,i),
						      at_or_last(bleaching_factors,i));
	  auto name = "marine_particles_"+names[i];
	  add_profile(m,name);
	}
      }
    }
    void add_marine_cdom() {
      std::vector<std::string> names = c_.get_vector<std::string>("mcdom_names");
      std::vector<double> scaling_factors = c_.get_vector<double>("mcdom_scaling_factors");
      ensure(names.size()==scaling_factors.size(), "marine cdom");
      for (size_t i = 0; i<names.size(); i++) {
	if (scaling_factors.at(i) > 0) {
	  auto m = std::make_shared<marine_cdom>(names.at(i), scaling_factors.at(i));
	  auto name = "marine_cdom_"+names[i];
	  add_profile(m,name);
	}
      }
    }
  private:
    template<class T>
    double at_or_last(const std::vector<T>& v, size_t i) {
      if (i > v.size()-1)
	return v.back();
      return v[i];
    }
    void ensure(bool b, const std::string& s) {
      if (not b)
	throw std::runtime_error("ocean "+s+" error");
    }
  };
}
}

#endif
