#ifndef flick_material_henyey_greenstein
#define flick_material_henyey_greenstein

#include "material.hpp"
#include "../polarization/rayleigh_mueller.hpp"

namespace flick {
namespace material {
  class henyey_greenstein : public basic_material {
    flick::absorption_coefficient ac_;
    flick::scattering_coefficient sc_;
    asymmetry_factor g_;
  public:
    henyey_greenstein(const flick::absorption_coefficient& ac,
		      const flick::scattering_coefficient& sc,
		      const flick::asymmetry_factor& g)
      : ac_{ac}, sc_{sc}, g_{g} {}
    double absorption_coefficient() {
      return ac_();
    }
    double scattering_coefficient() {
      return sc_();
    }
    auto mueller_matrix() {
      mueller m;
      m.add(0,0,flick::henyey_greenstein(g_()).phase_function(direction_.theta()));
      return m;
    }
    double refractive_index() {
      return 1;
    }
    friend std::ostream& operator<<(std::ostream &os,
				    henyey_greenstein& hg) {
      stream_basic_material(os,hg);
      os << ", asymmetry factor " << hg.g_();
      return os;
    }
  };
}
}

#endif
