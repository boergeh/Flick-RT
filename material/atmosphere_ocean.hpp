#ifndef flick_atmosphere_ocean
#define flick_atmosphere_ocean

#include "../material/ocean.hpp"
#include "../material/atmosphere.hpp"

namespace flick {
namespace material {
  struct atmosphere_ocean : public mixture {
    struct configuration : basic_configuration {
      configuration() {
	add_configuration(ocean::configuration());
	add_configuration(atmosphere::configuration());
      }
    };
    atmosphere_ocean(const basic_configuration& c=atmosphere_ocean::configuration())
      : mixture(range(0,constants::pi,c.get<size_t>("angles")).linspace(),
		atmospheric_state(c.get<size_t>("heights")).height_grid()) {
      should_update_iops(false);
      //add_material<ocean>(c);
      add_material<atmosphere>(c);
      should_update_iops(true);
      update_iops();
    }
  };
}
}

#endif
