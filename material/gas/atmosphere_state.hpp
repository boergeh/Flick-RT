#ifndef flick_atmosphere_state
#define flick_atmosphere_state

#include <set>
#include <map>
#include "profile.hpp"

namespace flick {
  class atmosphere_state {
    size_t n_points_;
    std::set<std::string> gases_{"co2","h2o","o2","o3"};
    std::map<std::string, concentration_profile> gas_concentrations_;
    std::map<std::string, double> gas_scaling_factors_;
    double temperature_scaling_factor_;
    double pressure_scaling_factor_;
    double air_scaling_factor_;
    concentration_profile air_concentration_{"air.txt",n_points_};
    temperature_profile temperature_{"temperature.txt",n_points_};
    pressure_profile pressure_{"pressure.txt",n_points_};
  public:
    atmosphere_state()
      : atmosphere_state(constants::T_stp,constants::P_stp) {};
    atmosphere_state(double surface_temperature, double surface_pressure,
		     size_t n_points=50)
      : n_points_{n_points} {
      load_gases();
      temperature_scaling_factor_ = surface_temperature / temperature_.value(0);
      pressure_scaling_factor_ = surface_pressure / pressure_.value(0);
      air_scaling_factor_ = pressure_scaling_factor_;
      for(auto g : gases_) {
	gas_scaling_factors_[g] = 1;
      }
    }
    std::vector<double> height_grid() const {
      return air_concentration_.grid();
    }
    void remove_gas(std::string gas_name) {
      gases_.erase(gas_name);
    }
    void scale_concentration(const std::string& gas_name, double factor)
    // Keeping total pressure
    {
      gas_scaling_factors_.at(gas_name) *= factor; 
    }
    void list_gases() const {
      for(auto g : gases_) {
	std::cout << g << "  "; 
      }
    }
    std::vector<std::string> gas_names() const {
      std::vector<std::string> s;
      for(auto g : gases_) {
	s.push_back(g);
      }
      return s;
    }
    double stp_thickness(const std::string& gas_name) const {
      return gas_concentrations_.at(gas_name).stp_thickness() *
	gas_scaling_factors_.at(gas_name);
    }
    double per_area(const std::string& gas_name) const {
      return gas_concentrations_.at(gas_name).per_area() *
	gas_scaling_factors_.at(gas_name);
    }
    double temperature(double height) const {
      return temperature_.value(height) * temperature_scaling_factor_;
    }
    double pressure(double height) const {
      return pressure_.value(height) * pressure_scaling_factor_;
    }
    double air_concentration(double height) const {
      return air_concentration_.value(height) * air_scaling_factor_;
    }
    double gas_concentration(const std::string& gas_name, double height) const {
      return gas_concentrations_.at(gas_name).value(height) *
	gas_scaling_factors_.at(gas_name);
    }
  private:
    void load_gases() {
      //air_concentration_ = concentration_profile("air.txt",n_points_);
      for(auto g : gases_) {
	gas_concentrations_[g] = concentration_profile(g+".txt",n_points_);
      }
    }
  };
}

#endif
  
