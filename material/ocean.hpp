#ifndef flick_material_ocean
#define flick_material_ocean

#include "water/pure_water.hpp"
#include "water/cdom.hpp"
#include "mixture.hpp"
#include "../environment/configuration.hpp"

namespace flick {
namespace material {
  struct ocean : public mixture {
    struct configuration : basic_configuration {
      configuration() {
	add<size_t>("angles", 20,
		    "Number of computational phase function angles");
	add<double>("bottom_depth", 200, "Depth of water column [m]");
	add<double>("cdom_440", 0.01, "CDOM absorption coefficient at 440 nm [1/m]");
	add<double>("cdom_slope", 0.017, "CDOM absorption spectrum slope [1/nm]");
	add<double>("pure_water_vf", 1, "Fraction of volume filled with pure water.");
      }
    };
  private:
    basic_configuration c_;
  public:
    ocean(const basic_configuration& c=ocean::configuration())
      : mixture(angle_range(c.get<size_t>("angles")),
		{-c.get<double>("bottom_depth"), 0}) {
      c_ = c;
      //should_update_iops(false);
      add_pure_water();
      add_cdom();
      //should_update_iops(true);
      //update_iops();
    }
  private:
    void add_pure_water() {
      add_material<pure_water>(30,290,c_.get<double>("pure_water_vf"));
    }
    void add_cdom() {
      add_material<cdom>(c_.get<double>("cdom_440"),c_.get<double>("cdom_slope"));
    }
  };
}
}

#endif
