#ifndef flick_atmosphere_state
#define flick_atmosphere_state

#include "profile.hpp"

namespace flick {
  enum gas {air,co2,h2o,no2,o2,o3};
  class atmosphere_state {
    const int n_gases_{6};
    std::vector<concentration_profile> gas_concentrations_;
    std::vector<double> gas_scaling_factors_;
    double temperature_scaling_factor_;
    double pressure_scaling_factor_;
    temperature_profile temperature_{"temperature.txt"};
    pressure_profile pressure_{"pressure.txt"};
  public:
    atmosphere_state()
      : atmosphere_state(constants::T_stp,constants::P_stp) {};
    atmosphere_state(double surface_temperature, double surface_pressure) {
      read_gases();
      temperature_scaling_factor_ = surface_temperature / temperature_.value(0);
      pressure_scaling_factor_ = surface_pressure / pressure_.value(0);
      gas_scaling_factors_ = std::vector<double>(n_gases_,pressure_scaling_factor_);
    }
    void scale_concentration(gas g, double factor)
    // Keeping total pressure
    {
      gas_scaling_factors_[g] *= factor; 
    }
    double stp_thickness(gas g) const {
      return gas_scaling_factors_[g] * gas_concentrations_[g].stp_thickness();
    }
    double temperature(double height) const {
      return temperature_.value(height) * temperature_scaling_factor_;
    }
    double pressure(double height) const {
      return pressure_.value(height) * pressure_scaling_factor_;
    }
  private:
    void read_gases() {
      gas_concentrations_.push_back(concentration_profile("air.txt"));
      gas_concentrations_.push_back(concentration_profile("co2.txt"));
      gas_concentrations_.push_back(concentration_profile("h2o.txt"));
      gas_concentrations_.push_back(concentration_profile("no2.txt"));
      gas_concentrations_.push_back(concentration_profile("o2.txt"));
      gas_concentrations_.push_back(concentration_profile("o3.txt"));
    }
  };
}

#endif
  
