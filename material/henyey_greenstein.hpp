#ifndef flick_material_henyey_greenstein
#define flick_material_henyey_greenstein

#include "monocrome_iop.hpp"

namespace flick {
namespace material {
  class henyey_greenstein : public monocrome_iop {
    flick::asymmetry_factor g_;
  public:
    henyey_greenstein(const flick::absorption_coefficient& ac,
		      const flick::scattering_coefficient& sc,
		      const flick::asymmetry_factor& g)
      : monocrome_iop{ac,sc}, g_{g} {
      set(flick::asymmetry_factor{g_()});
    }  
    void refractive_index(double ri) {
      refractive_index_ = ri;
    }
    double absorption_coefficient() {
      return ac_();
    }
    double scattering_coefficient() {
      return sc_();
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) {
      mueller m;
      double theta = angle(scattering_direction);
      m.add(0,0,flick::henyey_greenstein(g_()).phase_function(theta));
      return m;
    }
    double refractive_index() {
      return refractive_index_;;
    }
  };
}
}

#endif
