#ifndef flick_material_henyey_greenstein
#define flick_material_henyey_greenstein

#include "material.hpp"
#include "../polarization/rayleigh_mueller.hpp"
#include "named_bounded_types.hpp"

namespace flick {
namespace material {
  class henyey_greenstein : public base {
    flick::absorption_coefficient ac_;
    flick::scattering_coefficient sc_;
    flick::asymmetry_factor g_;
  public:
    henyey_greenstein(const flick::absorption_coefficient& ac,
		      const flick::scattering_coefficient& sc,
		      const flick::asymmetry_factor& g)
      : ac_{ac}, sc_{sc}, g_{g} {}
    
    henyey_greenstein(double ac,
		      double sc,
		      double g)
      : ac_{ac}, sc_{sc}, g_{g} {}
    
    double absorption_optical_depth(double distance) {
      return absorption_coefficient()*distance;
    }
    double scattering_optical_depth(double distance) {
      return scattering_coefficient()*distance;
    }
    double absorption_distance(double absorption_optical_depth) {
       double ac = absorption_coefficient();
      if (ac > 0)
	return absorption_optical_depth/ac;
      return std::numeric_limits<double>::max(); 
    }
    double scattering_distance(double scattering_optical_depth) {
      double sc = scattering_coefficient();
      if (sc > 0)
	return scattering_optical_depth/sc;
      return std::numeric_limits<double>::max(); 
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
      return 1;
    }
    double sampling_asymmetry_factor() {
      return g_();
    }
    friend std::ostream& operator<<(std::ostream &os,
				    henyey_greenstein& hg) {
      //stream_basic_material(os,hg);
      os << ", asymmetry factor " << hg.g_();
      return os;
    }
  };
}
}

#endif