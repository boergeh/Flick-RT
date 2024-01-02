#ifndef flick_sun_position
#define flick_sun_position

#include <numbers>
#include "time_point.hpp"

namespace flick {
  class earth_orbit
  // https://en.wikipedia.org/wiki/Position_of_the_Sun
  // https://en.wikipedia.org/wiki/Equation_of_time
  {
    int year_;
    double day_of_year_;
    double to_radians_ = std::numbers::pi/180;
  public:
    earth_orbit(int year, double day_of_year)
      : year_{year}, day_of_year_{day_of_year} {
    }
    double distance() const {
      const double au = 149597870700;
      double g = (357.528 + 0.9856003*days_since_2000())*to_radians_;
      return (1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g))*au;
    }
    double equation_of_time() const {
      double D = 6.24004077+0.01720197*(365.25*(year_-2000)+day_of_year_);
      double delta_t =  -7.659*sin(D)+9.863*sin(2*D+3.5932); 
      return delta_t * 60; // [seconds]
    }
    double declination() const {
      double N = days_since_2000();
      double u = 0.98565*(N-2)*to_radians_;
      double v = (0.98565*(N+10)+1.914*sin(u))*to_radians_;
      return -asin(0.39779*cos(v));
    }
  private:
    double days_since_2000() const {
      return julian_date_since_2000(time_point(year_,1,1,0,0,0))+day_of_year_;
    }
  };
  
  class sun_position
  // https://en.wikipedia.org/wiki/Solar_azimuth_angle
  {
    time_point time_point_; // UTC
    double latitude_; // [radians]
    double longitude_; // [radians]
    earth_orbit earth_orbit_;
    std::vector<double> S_;
    double to_radians_ = std::numbers::pi/180;
  public:
    sun_position(const time_point& t, double latitude, double longitude)
      : time_point_{t}, latitude_{latitude}, longitude_{longitude},
       	earth_orbit_{t.year(),day_of_year(t)} {
      S_ = sun_direction();      
    }
    double zenith_angle() const {
      return acos(S_[2]);
    }
    double azimuth_angle() const {
      return std::atan2(-S_[0],-S_[1]);
    }
  private:
    double latitude_of_subsolar_point() const {
      return earth_orbit_.declination();
    }
    double longitude_of_subsolar_point() const {
      double E_min = earth_orbit_.equation_of_time()/60/60; //[hours]
      double T_utc = hour_of_day(time_point_);
      double degrees_per_hour = 15;
      return -degrees_per_hour*(T_utc-12+E_min)*to_radians_;
    }
    std::vector<double> sun_direction() const {
      double phi_s = latitude_of_subsolar_point();
      double phi_0 = latitude_;
      double lambda_s = longitude_of_subsolar_point();
      double lambda_0 = longitude_;
      double dl = lambda_s-lambda_0;
      std::vector<double> S(3);
      S[0] = cos(phi_s)*sin(dl);
      S[1] = cos(phi_0)*sin(phi_s)-sin(phi_0)*cos(phi_s)*cos(dl);
      S[2] = sin(phi_0)*sin(phi_s)+cos(phi_0)*cos(phi_s)*cos(dl);
      return S;
    }
  };
}

#endif
