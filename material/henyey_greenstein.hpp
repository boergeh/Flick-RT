#ifndef flick_material_henyey_greenstein
#define flick_material_henyey_greenstein

#include "monocrome_iop.hpp"

namespace flick {
namespace material {
  class henyey_greenstein : public monocrome_iop {
  public:
    using monocrome_iop::monocrome_iop;
    mueller mueller_matrix(const unit_vector& scattering_direction) {
      mueller m;
      double theta = angle(scattering_direction);
      m.add(0,0,flick::henyey_greenstein(g_()).phase_function(theta));
      return m;
    }
  };
}
}

#endif
