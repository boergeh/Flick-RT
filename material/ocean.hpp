#ifndef flick_material_ocean
#define flick_material_ocean

#include "water/pure_water.hpp"
#include "mixture.hpp"
#include "../environment/configuration.hpp"

namespace flick {
namespace material {
  struct ocean : public mixture {
    struct configuration : basic_configuration {
      configuration() {
	//add<size_t>("angles", 2,
	//	    "Number of computed phase function angles");
	add<double>("bottom_depth", 200, "[m]");
      }
    };

    ocean(const basic_configuration& c=ocean::configuration())
      : mixture(range(0,constants::pi,50).linspace(),
		       {-c.get<double>("bottom_depth"), 0}) {
      should_update_iops(false);
      add_pure_water();     
      should_update_iops(true);
      update_iops();
    }
  private:
    void add_pure_water() {
      add_material<pure_water>();
    }
  };
}
}

#endif
