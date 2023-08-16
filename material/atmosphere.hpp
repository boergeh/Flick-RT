#ifndef flick_material_atmosphere
#define flick_material_atmosphere

#include "gas/air.hpp"
#include "gas/atmospheric_state.hpp"
#include "aerosols/aerosols.hpp"
#include "mixture.hpp"

namespace flick {
namespace material {
  struct atmosphere : public mixture {
    atmosphere(size_t n_angles, const std::vector<double>& heights)
      : mixture(range(0,constants::pi,n_angles).linspace(), heights) {     
      set_material<air>(atmospheric_state());
      //set_material<rural_aerosols>(0.1, 0.5);
      update_iops();
    }    
  };
}
}

#endif
