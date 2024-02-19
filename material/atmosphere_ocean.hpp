#ifndef flick_atmosphere_ocean
#define flick_atmosphere_ocean

#include "../material/ocean.hpp"
#include "../material/atmosphere.hpp"

namespace flick {
namespace material {
  struct atmosphere_ocean : public mixture<pe_function> {
    struct configuration : basic_configuration {
      configuration() {
	add_configuration(ocean::configuration());
	add_configuration(atmosphere::configuration());
      }
    };
    atmosphere_ocean(const basic_configuration& c=atmosphere_ocean::configuration())
      : mixture(angle_range(c.get<size_t>("n_angles")), height_grid(c)) {
      auto_update_iops(false);
      add_material<ocean>(c);
      size_t n_oce = ocean::height_grid(c).size();
      size_t n_atm = atmosphere::height_grid(c).size();      
      set_range<ocean>(0, n_oce-1);
      add_material<atmosphere>(c);
      set_range<atmosphere>(n_oce, n_oce+n_atm-1);
      auto_update_iops(true);
    }
    static stdvector height_grid(const basic_configuration& c) {
      stdvector a = ocean::height_grid(c);
      stdvector b = atmosphere::height_grid(c);
      a.insert(a.end(),b.begin(),b.end());
      return a;
    }
  };
}
}

#endif
