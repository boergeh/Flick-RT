#ifndef flick_profile
#define flick_profile

#include <fstream>
#include <set>
#include "../../numeric/constants.hpp"
#include "../../numeric/function.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
  class profile {
  protected:
    size_t n_points_;
    double column_integral_;
    pe_function f_;
  public:
    profile() = default;
    profile(const std::string& file_name, size_t n_points = 50)
      : n_points_{n_points} {
      f_ = read<pe_function>
	("material/gas/profile_input/"+file_name);
      column_integral_ = f_.integral(0,120e3);
      f_ = sparse_profile_distribution();
    }
    double value(double height_above_surface) const {
      return f_.value(height_above_surface)*column_integral_;
    }    
    std::vector<double> grid() const {
      return f_.x();
    }   
  private:
     pe_function sparse_profile_distribution() const {
      std::set<size_t> sg = sparse_grid();
      pe_function f;
      auto x = f_.x();
      auto y = f_.y();
      for (auto i : sg) {
	f.append({x[i],y[i]});
      }
      return f.scale_y(1/f.integral());
    }
    std::set<size_t> sparse_grid() const {
      std::vector<double> gd =
	range(0,f_.size()-1,n_points_).linspace();
      std::set<size_t> gi;
      for (auto i : gd) {
	gi.insert(i);
      }
      return gi;
    }   
  };
  
  struct concentration_profile : public profile {
    using profile::profile;
    double per_area() const {
      return column_integral_;
    }
    double stp_thickness() const {
      return per_area() * constants::k_B * constants::T_stp
	/ constants::P_stp;
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
  
