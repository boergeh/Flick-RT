#ifndef flick_time_point
#define flick_time_point

#include <ostream>

namespace flick {
  class time_point
  // See https://en.wikipedia.org/wiki/Julian_day
  {
    int year_;
    int month_;
    int day_;
    int hour_;
    int minute_;
    double second_;
  public:
    time_point() : time_point(-4713,11,24,12,0,0) {
    }
    time_point(int year, int month, int day,
	       int hour, int minute, double second)
      : year_{year}, month_{month}, day_{day},
	hour_{hour}, minute_{minute}, second_{second} {}
    double julian_date() const {
      return julian_day_number()+(hour_-12)/24.+minute_/1440.+second_/86400.;
    }
    int year() const {
      return year_;
    }
    int month() const {
      return month_;
    }
    int day() const {
      return day_;
    }
  private:
    int julian_day_number() const {
      int Y = year_;
      int M = month_;
      int D = day_;
      return (1461 * (Y + 4800 + (M - 14)/12))/4 +
	(367 * (M - 2 - 12 * ((M - 14)/12)))/12 -
	(3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075;
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const time_point& t) {
      os << t.year_ << " " << t.month_ << " " << t.day_ << " "
	 << t.hour_ << " " << t.minute_ << " " << t.second_;
      return os;
    }
  };
  
  double day_of_year(const time_point& t) {
    return t.julian_date() -
      time_point(t.year(),1,1,0,0,0).julian_date();
  }
  double hour_of_day(const time_point& t) {
    double days = t.julian_date() -
      time_point(t.year(),t.month(),t.day(),0,0,0).julian_date();
    return days*24;
  }
  double julian_date_since_2000(const time_point& t) {
    time_point J2000 = {2000,1,1,12,0,0};
    return t.julian_date() - J2000.julian_date();
  }
}

#endif
