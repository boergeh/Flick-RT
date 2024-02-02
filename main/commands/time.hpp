#include "basic_command.hpp"
#include "../../astronomy/time_point.hpp"

namespace flick {
  namespace command {
    struct time : public basic_command {
      time():basic_command("time"){};
      void run() {
	std::string epoch = a(1);
	if (size()==3) {
	  std::cout << std::fixed<<std::setprecision(0);
	  double days_TT = std::stod(a(2));
	  if (epoch=="julian_date")
	    std::cout << time_converter().TT_to_UTC(make_time_point(days_TT));
	  else if (epoch=="datenum")
	    std::cout <<  time_converter().TT_to_UTC(make_time_point(datenum_to_julian_date(days_TT)));
	  else if (epoch=="J2000")
	    std::cout <<  time_converter().TT_to_UTC(make_time_point(J2000_to_julian_date(days_TT)));
	  else
	    error();
	} else if (size()==8) {
	  std::cout << std::fixed<<std::setprecision(5);
	  flick::time_point tp_utc = get_time_point();
	  flick::time_point t = time_converter().UTC_to_TT(tp_utc); 
	  if (epoch=="julian_date")
	    std::cout << t.julian_date();
	  else if (epoch=="datenum")
	    std::cout << t.matlab_datenum();
	  else if (epoch=="J2000")
	    std::cout << t.J2000();
	  else
	    error();
	} else {
	  error();
	}
	std::cout << std::defaultfloat<<std::setprecision(6);
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
