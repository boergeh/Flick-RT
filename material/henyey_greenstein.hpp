#ifndef flick_material_henyey_greenstein
#define flick_material_henyey_greenstein

#include "monocrome_iop.hpp"

namespace flick {
namespace material {
  class henyey_greenstein : public monocrome_iop {
  public:
    using monocrome_iop::monocrome_iop;
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double theta = angle(scattering_direction);
      m.add(0,0,flick::henyey_greenstein(g_()).phase_function(theta));
      return m;
    }
  };
  class white_isotropic : public henyey_greenstein {
  public:
    white_isotropic(double scat_coef)
      : henyey_greenstein(flick::absorption_coefficient{0},
			  flick::scattering_coefficient{scat_coef},
			  flick::asymmetry_factor{0}) {}
  };
}
}

#endif
