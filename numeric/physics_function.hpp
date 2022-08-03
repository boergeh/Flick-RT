#ifndef flick_physics_functions
#define flick_physics_functions

#include "constants.hpp"
#include "function.hpp"
#include "range.hpp"
#include "../environment/exception.hpp"

namespace flick {
  class henyey_greenstein
  // Henyey-Greenstein scattering phase function
  {
    double asymmetry_factor_{0};
  public:
    henyey_greenstein(double g) : asymmetry_factor_{g} {}
    double phase_function(double angle)
    // Note that integral over 4*pi equals one 
    {
      double g = asymmetry_factor_;
      double pi = constants::pi;
      double arg = 1+pow(g,2)-2*g*cos(angle);
      return 1/(4*pi)*(1-pow(g,2))/pow(arg,3./2);
    }
    double inverted_accumulated_angle(double fraction)
    // Note that fraction has range [0 1]      
    {
      double g = asymmetry_factor_;
      if (fabs(g) < 1e-9) // isotropic
	return acos(1-2*fraction);
      double arg = (1-pow(g,2))/(1-g+2*g*fraction);
      double mu = (1+pow(g,2)-pow(arg,2))/(2*g);
      return acos(-mu);
    }
  };
}

#endif
