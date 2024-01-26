#ifndef flick_time_point
#define flick_time_point

#include <iostream>

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
    struct epoch {
      static time_point matlab_datenum() {
	return time_point(-1,12,31,0,0,0);
      }
      static time_point julian_day() {
	return time_point(-4713,11,24,12,0,0);
      }
      static time_point J2000() {
	return time_point(2000,1,1,12,0,0);
      }
    };
    time_point() : time_point(epoch::julian_day()) {
    }
    time_point(int year, int month, int day,
	       int hour, int minute, double second)
      : year_{year}, month_{month}, day_{day},
	hour_{hour}, minute_{minute}, second_{second} {
    }
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
    int hour() const {
      return hour_;
    }
    int minute() const {
      return minute_;
    }
    double second() const {
      return second_;
    }
    double day_of_year() const {
      return julian_date() - time_point(year_,1,1,0,0,0).julian_date();
    }
    double hour_of_day() const {
      double days = julian_date() -
	time_point(year_,month_,day_,0,0,0).julian_date();
      return days*24;
    }
    double J2000() const {
      return julian_date() - epoch::J2000().julian_date();
    }
    double matlab_datenum() const {
      return julian_date() - epoch::matlab_datenum().julian_date();
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

  double datenum_to_julian_date(double matlab_datenum) {
    return matlab_datenum + time_point::epoch::matlab_datenum().julian_date();
  }
  
  double J2000_to_julian_date(double J2000) {
    return J2000 + time_point::epoch::J2000().julian_date();
  }

  time_point make_time_point(double julian_date)
  // https://en.wikipedia.org/wiki/Julian_day
  {
    int y = 4716;
    int j = 1401;
    int m = 2;
    int n = 12;
    int r = 4;
    int p = 1461;
    int v = 3;
    int u = 5;
    int s = 153;
    int w = 2;
    int B = 274277;
    int C = -38;      
    int J = round(julian_date);
    int f = J+j+(((4*J+B)/146097)*3)/4+C;
    int e = r*f+v;
    int g = (e % p)/r;
    int h = u*g+w;
    int D = (h % s)/u+1;
    int M = (h/s+m) % n + 1;
    int Y = e/p-y+(n+m-M)/n;
    int hr = (julian_date-time_point(Y,M,D,0,0,0).julian_date())*24;
    int mi = (julian_date-time_point(Y,M,D,hr,0,0).julian_date())*60*24;
    double se = (julian_date-time_point(Y,M,D,hr,mi,0).julian_date())*60*60*24;
    return time_point(Y,M,D,hr,mi,se);
  }
}

#endif
