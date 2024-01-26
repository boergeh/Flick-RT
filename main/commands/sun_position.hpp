#include "basic_command.hpp"
#include "../../astronomy/sun_position.hpp"

namespace flick {
  namespace command {
    struct sun_position : public basic_command {
      double to_radians = std::numbers::pi/180;
      double to_degrees = 180/std::numbers::pi;
      sun_position():basic_command("sun_position"){};
      void run() {
	if (a(1)=="distance") {
	  flick::time_point t = get_time_point(); 
	  std::cout << flick::earth_orbit(t.year(),t.day_of_year()).distance()/constants::au;
	}
	else if (a(1)=="zenith_angle") {
	  double latitude = std::stod(a(8))*to_radians;
	  double longitude = std::stod(a(9))*to_radians;
	  std::cout << flick::sun_position(get_time_point(),latitude,longitude).zenith_angle()*to_degrees;
	}
	else if (a(1)=="azimuth_angle") {
	  double latitude = std::stod(a(8))*to_radians;
	  double longitude = std::stod(a(9))*to_radians;
	  std::cout << flick::sun_position(get_time_point(),latitude,longitude).azimuth_angle()*to_degrees;
	}
	else {
	  error();
	}
      }
      time_point get_time_point() {
	int year = std::stoi(a(2));
	int month = std::stoi(a(3));
	int day = std::stoi(a(4));
	int hour = std::stoi(a(5));
	int minute = std::stoi(a(6));
	double second = std::stod(a(7));
	return time_point(year,month,day,hour,minute,second);
      }
    };
  }
}
