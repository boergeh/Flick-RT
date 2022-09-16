#ifndef flick_material_monocrome_iop
#define flick_material_monocrome_iop

#include "material.hpp"

namespace flick {
namespace material {
  class monocrome_iop : public base {
  protected:
    flick::absorption_coefficient ac_;
    flick::scattering_coefficient sc_;
    flick::asymmetry_factor g_;
    double real_refractive_index_;
  public:
    monocrome_iop(const flick::absorption_coefficient& ac,
		  const flick::scattering_coefficient& sc,
		  const flick::asymmetry_factor& g,
		  double real_refractive_index = 1)
      : ac_{ac}, sc_{sc}, g_{g}, real_refractive_index_{real_refractive_index} {}    
    double absorption_coefficient() {
      return ac_();
    }
    double scattering_coefficient() {
      return sc_();
    }
    double asymmetry_factor() {
      return g_();
    }
    double real_refractive_index() {
      return real_refractive_index_;
    }
  };
  
}
}

#endif
