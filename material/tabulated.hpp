#ifndef flick_material_tabulated
#define flick_material_tabulated

#include "monocrome_iop.hpp"
#include "../numeric/function.hpp"

namespace flick {
namespace material {
  class tabulated : public monocrome_iop {
    tabulated_phase_function p_;
  public:
    tabulated(const flick::absorption_coefficient& ac,
	      const flick::scattering_coefficient& sc,
	      const tabulated_phase_function& p,
	      double real_refractive_index = 1)
      : monocrome_iop{ac,sc,flick::asymmetry_factor{p.asymmetry_factor()},
      real_refractive_index}, p_{p} {
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double theta = angle(scattering_direction);
      m.add(0,0,p_.value(theta));
      return m;
    }
  };
}
}

#endif
