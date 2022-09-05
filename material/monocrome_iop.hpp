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
    //std::complex<double> refractive_index_;
    double refractive_index_{1};
  public:
    monocrome_iop(const flick::absorption_coefficient& ac,
		  const flick::scattering_coefficient& sc)
      : ac_{ac}, sc_{sc} {}    
    void refractive_index(double ri) {
      refractive_index_ = ri;
    }
    double absorption_coefficient() {
      return ac_();
    }
    double scattering_coefficient() {
      return sc_();
    }
    double refractive_index() {
      return refractive_index_;;
    }
  };
}
}

#endif
