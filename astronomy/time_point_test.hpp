#include "time_point.hpp"
#include "../environment/unit_test.hpp"

namespace flick {
  begin_test_case(time_point_test_A) {
    time_point J0 = {-4713,11,24,12,0,0};
    check_small(J0.julian_date());
    check_small(time_point().julian_date());
    time_point  tp = {2000,1,1,12,0,0};
    check_small(tp.julian_date()-2451545.0);
    check_small(tp.J2000());
  } end_test_case()
  
  begin_test_case(time_point_test_B) {
    time_point tp = {2024,1,25,13,37,0};
    check_small(tp.matlab_datenum()-739276.5673611111,1e-9);
  } end_test_case()
  
  begin_test_case(time_point_test_C) {
    time_point tp1 = {2024,1,25,13,37,1};
    double jd = tp1.julian_date();
    time_point tp2 = make_time_point(jd);
    std::cout << tp2;
    check(tp1.year()==tp2.year());
    check(tp1.month()==tp2.month());
    check(tp1.day()==tp2.day());
    check(tp1.hour()==tp2.hour());
    check(tp1.minute()==tp2.minute());
    check_close(tp1.second(),tp2.second(),0.0006_pct);
  } end_test_case()
  
  begin_test_case(time_point_test_D) {
    time_point tp1 = {2024,1,25,13,37,1};
    double dn = tp1.matlab_datenum();
    time_point tp2 = make_time_point(datenum_to_julian_date(dn));
    check(tp1.year()==tp2.year());
    check(tp1.month()==tp2.month());
    check(tp1.day()==tp2.day());
    check(tp1.hour()==tp2.hour());
    check(tp1.minute()==tp2.minute());
    check_close(tp1.second(),tp2.second(),0.0006_pct);
  } end_test_case()
  
  begin_test_case(time_point_test_E) {
    time_point tp1 = {2024,1,26,0,0,1};
    time_point tp2 = make_time_point(tp1.julian_date());
    std::cout << tp2;
    check(tp1.year()==tp2.year());
    check(tp1.month()==tp2.month());
    check(tp1.day()==tp2.day());
    check(tp1.hour()==tp2.hour());
    check(tp1.minute()==tp2.minute());
    check_close(tp1.second(),tp2.second(),0.0006_pct);
  } end_test_case()
}
