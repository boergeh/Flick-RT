#ifndef flick_material_ocean
#define flick_material_ocean

#include "water/pure_water.hpp"
#include "water/cdom.hpp"
#include "water/phytoplankton.hpp"
#include "water/nap.hpp"
#include "marine_particles/marine_particles.hpp"
#include "mixture.hpp"
#include "../environment/configuration.hpp"

namespace flick {
namespace material {
  struct ocean : public mixture {
    struct configuration : basic_configuration {
      configuration() {
	add<size_t>("n_angles", 20, "Number of computational phase function angles");
	add<double>("bottom_depth", 200, "Depth of water column [m]");
	add<double>("cdom_440", 0.01, "CDOM absorption coefficient at 440 nm [1/m]");
	add<double>("cdom_slope", 0.017, "CDOM absorption spectrum slope [1/nm]");
	add<double>("chl_concentration", 1e-6, "Chlorophyll concentration in the water column [kg/m^3]");
	add<double>("nap_concentration", 1e-3, "Nonalgal particle dry mass concentration in the water column [kg/m^3]");
	add<std::string>("mp_names", "MP21_PA61", "A list of names of additional tabulated marine particles");
	add<double>("mp_concentrations", 0, R"(A list of dry mass concentrations of additional tabulated marine
particles in the water column [kg/m^3], one concentration value for
each name)");
      }
    };
  private:
    basic_configuration c_;
  public:
    ocean(const basic_configuration& c=ocean::configuration())
      : mixture(angle_range(c.get<size_t>("n_angles")),
		{-c.get<double>("bottom_depth"), 0}) {
      c_ = c;
      add_pure_water();
      add_cdom();
      add_phytoplankton();
      add_nap();
      add_marine_particles();
    }
  private:
    void add_pure_water() {
      add_material<pure_water>(30,290);
    }
    void add_cdom() {
      add_material<cdom>(c_.get<double>("cdom_440"),c_.get<double>("cdom_slope"));
    }
    void add_phytoplankton() {
      add_material<phytoplankton>(c_.get<double>("chl_concentration"));
    }
    void add_nap() {
      add_material<nap>(c_.get<double>("nap_concentration"));
    }
    void add_marine_particles() {
      std::vector<std::string> names = c_.get_vector<std::string>("mp_names");
      std::vector<double> concentrations = c_.get_vector<double>("mp_concentrations");
      for (size_t i = 0; i<names.size(); i++) {
	name_extension(std::to_string(i));
	add_material<marine_particles>(names.at(i), concentrations.at(i));
	name_extension("");
      }
    }
  };
}
}

#endif
