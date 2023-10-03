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
      : mixture(angle_range(c.get<size_t>("angles")), boundaries()) {
      should_update_iops(false);
      add_material<ocean>(c);
      set_range<ocean>(0,1);
      add_material<atmosphere>(c);
      set_range<atmosphere>(2,boundaries().size()-1);
      should_update_iops(true);
      update_iops();
    }
  private:
    stdvector boundaries() const {
      stdvector h1 = {-200, -0.001};
      //stdvector h2 = atmospheric_state(c.get<size_t>("heights")).height_grid();
      stdvector h2 = atmospheric_state(8).height_grid();
      h1.insert(h1.end(),h2.begin(),h2.end());
      return h1;
    }
  };
}
}

#endif
