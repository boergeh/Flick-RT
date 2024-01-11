#ifndef flick_constants
#define flick_constants

namespace flick {
  namespace constants
  {
    inline constexpr double pi = 3.141592653589793238462643383279;
    
    // Boltzman constant (exact)
    inline constexpr double k_B = 1.380649e-23; // [J/K]

    // Planck constant (exact)
    inline constexpr double h = 6.62607015e-34; // [J s]

    // Speed of light in vacuum (exact)
    inline constexpr double c = 299792458; // [m/s]

    // Astronomical unit (exact)
    inline constexpr double au = 149597870700; // [m]

    // Standard atmospheric pressure
    inline constexpr double P_stp = 1e5; // [Pa]

    // Standard atmospheric temperature
    inline constexpr double T_stp = 273.15; // [K]

    // Natural atmospheric pressure
    inline constexpr double P_ntp = 101.325e3; // [Pa]

    // Natural atmospheric temperature
    inline constexpr double T_ntp = 293.15; // [K]

    // Degree Celcius for a absolute zero
    inline constexpr double zero_K = -273.15; // [Celsius]

    // Stefan-Boltzmann constant
    inline constexpr double sigma = 5.670374419184429453970996731889e-8; // [W/m^2/K^4]

    double to_celsius(double kelvin) {
      return kelvin + zero_K;
    }
    double to_kelvin(double celsius) {
      return celsius - zero_K;
    }
  }
}

#endif
