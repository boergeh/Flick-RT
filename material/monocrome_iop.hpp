#ifndef flick_material_monocrome_iop
#define flick_material_monocrome_iop

#include "material.hpp"
#include "../numeric/named_bounded_types.hpp"

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
    monocrome_iop(double ac, double sc, double g,
		  double real_refractive_index = 1)
      : ac_{flick::absorption_coefficient(ac)}, sc_{flick::scattering_coefficient{sc}},
	g_{flick::asymmetry_factor{g}}, real_refractive_index_{real_refractive_index} {}    
    double absorption_coefficient() const {
      return ac_();
    }
    double scattering_coefficient() const {
      return sc_();
    }
    double asymmetry_factor() const {
      return g_();
    }
    double real_refractive_index() const {
      return real_refractive_index_;
    }
  }; 
}
}

#endif
