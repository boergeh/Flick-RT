#include "basic_command.hpp"
#include "../../astronomy/time_point.hpp"

namespace flick {
  namespace command {
    struct time : public basic_command {
      time():basic_command("time"){};
      void run() {
	std::cout << std::setprecision(12);
	std::string epoch = a(1);
	if (size()==3) {
	  double days = std::stod(a(2));
	  if (epoch=="julian_date")
	    std::cout << make_time_point(days);
	  else if (epoch=="datenum")
	    std::cout << make_time_point(datenum_to_julian_date(days));
	  else
	    error();
	} else if (size()==8) { 
	  flick::time_point t = get_time_point();
	  if (epoch=="julian_date")
	    std::cout << t.julian_date();
	  else if (epoch=="datenum")
	    std::cout << t.matlab_datenum();
	  else
	    error();
	} else {
	  error();
	}
	std::cout << std::setprecision(6);
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
