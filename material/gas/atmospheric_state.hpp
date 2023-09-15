#ifndef flick_atmospheric_state
#define flick_atmospheric_state

#include <set>
#include <map>
#include "profile.hpp"

namespace flick {
  class atmospheric_state {
    size_t n_points_;
    std::set<std::string> gases_{"co2","h2o","o2","o3","no2"};
    std::map<std::string, concentration_profile> gas_concentrations_;
    std::map<std::string, double> gas_scaling_factors_;
    concentration_profile air_concentration_{"air.txt", n_points_};
    temperature_profile temperature_{"temperature.txt", n_points_};
    pressure_profile pressure_{"pressure.txt", n_points_};
    double temperature_scaling_factor_;
    double pressure_scaling_factor_;
    double air_scaling_factor_;
  public:
    atmospheric_state()
      : atmospheric_state(constants::T_stp, constants::P_stp, 50) {};
    atmospheric_state(size_t n_points)      
      : atmospheric_state(constants::T_stp, constants::P_stp, n_points) {};
    atmospheric_state(double surface_temperature, double surface_pressure,
		     size_t n_points=50)
      : n_points_{n_points} {
      temperature_scaling_factor_ = surface_temperature / temperature_.surface_value();
      pressure_scaling_factor_ = surface_pressure / pressure_.surface_value();
      air_scaling_factor_ = pressure_scaling_factor_;
      for(auto g : gases_) {
	gas_concentrations_[g] = concentration_profile(g+".txt", n_points_);
	gas_scaling_factors_[g] = air_scaling_factor_;
      }
    }
    std::vector<double> height_grid() const {
      return air_concentration_.grid();
    }
    void add_gas(const std::string& gas_name) {
      gases_.insert(gas_name);
    }
    void remove_all_gases() {
      gases_.clear();
    }
    void remove_gas(const std::string& gas_name) {
      gases_.erase(gas_name);
    }
    void scale_concentration(const std::string& gas_name, double factor)
    // Conserving total pressure
    {
      gas_scaling_factors_.at(gas_name) *= factor; 
    }
    void scale_to_stp_thickness(const std::string& gas_name, double new_thickness)
    {
      scale_concentration(gas_name, new_thickness / stp_thickness(gas_name));
    }
    void scale_to_fraction(const std::string& gas_name, double new_fraction)
    {
      scale_concentration(gas_name, new_fraction / fraction(gas_name));
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
    double fraction(const std::string& gas_name) const {
      return per_area(gas_name) / (air_concentration_.per_area() * air_scaling_factor_);
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
  };
}

#endif
  
