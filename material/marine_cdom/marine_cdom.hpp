#ifndef flick_material_marine_cdom
#define flick_material_marine_cdom

#include "../material.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class marine_cdom : public base {
    pl_function a_;
    double a440_;
  public:
    marine_cdom(const std::string& name, double scaling_factor) {
      std::string p = path()+"/material/marine_cdom/iop_tables";
      a_ = read<pl_function>(name+"_a.txt", p);

      // Subtract rayleigh attenuation
      /*
      std::vector<double> x = a_.x();
      std::vector<double> y(x.size());
      double x0 = 440;
      for (size_t i=0; i<x.size(); i++) {
	y[i] = a_.value(x[i]) - 0.5*a_.value(x0)*pow(x0/x[i],4);
      }
      */
      //a_ = pl_function(x,y);
      //a_.scale_y(0.8);
      
      a_.scale_y(scaling_factor);
      a_.add_constant_extrapolation();

      
    }
    double absorption_coefficient() const {
      const double to_nm = 1e9;
      return a_.value(wavelength()*to_nm);
    }
    double scattering_coefficient() const {
      return 0;
    }
  };
  
  template<int n>
  struct listable_marine_cdom : public marine_cdom {
    using marine_cdom::marine_cdom;
  };
}
}

#endif
