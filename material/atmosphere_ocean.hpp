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
      : mixture(angle_range(c.get<size_t>("angles")), boundaries(c.get<double>("bottom_depth"),c.get<size_t>("heights"))) {
      //should_update_iops(false);
      //c_ = c;
      add_material<ocean>(c);
      set_range<ocean>(0,1);
      add_material<atmosphere>(c);
      set_range<atmosphere>(2,boundaries(c.get<double>("bottom_depth"),c.get<size_t>("heights")).size()-1);
      //should_update_iops(true);
      update_iops();
    }
  private:
    //basic_configuration c_;
    stdvector boundaries(double bottom_depth, size_t n_atm_heights) const {
      stdvector h1 = {-bottom_depth, -1e-6};
      stdvector h2 = atmospheric_state(n_atm_heights).height_grid();
      h1.insert(h1.end(),h2.begin(),h2.end());
      return h1;
    }
  };
}
}

#endif