#ifndef flick_material_fournier_forand
#define flick_material_fournier_forand

#include "monocrome_iop.hpp"

namespace flick {
namespace material {
  class fournier_forand : public monocrome_iop {
  public:
    using monocrome_iop::monocrome_iop;
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double theta = angle(scattering_direction);
      m.add(0,0,flick::fournier_forand(g_()).value(cos(theta)));
      return m;
    }
  };
}
}

#endif
