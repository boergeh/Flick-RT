#ifndef flick_profile
#define flick_profile

#include <fstream>
#include "../../numeric/function.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
  class profile {
  protected:
    pe_function concentration_;
  public:
    profile(const std::string& file_name) {
      concentration_ = read<pe_function>
	("material/gas/profile_input/"+file_name);
    }
    double value(double height_above_surface) {
      return concentration_.value(height_above_surface);
    }
  };
  struct concentration_profile : public profile {
    using profile::profile;
    double stp_thickness() {
      double molecules_per_area = concentration_.integral(0,120e3);
      return molecules_per_area * constants::k_B * constants::T_stp / constants::P_stp;
    }
  };
  struct temperature_profile : public profile {
    using profile::profile;
  };
  struct pressure_profile : public profile {
    using profile::profile;
  };
}

#endif
  
