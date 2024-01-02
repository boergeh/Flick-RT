#include "sun_position.hpp"
#include "../environment/unit_test.hpp"

namespace flick {
  begin_test_case(sun_position_test_A) {
    double to_radians = std::numbers::pi/180;
    double to_degrees = 180/std::numbers::pi;

    time_point tp_utc = {2024,1,1,22,0,0};
    check_close(day_of_year(tp_utc),22./24,1e-7_pct);
    check_close(hour_of_day(tp_utc),22,1e-7_pct);
    
    earth_orbit eo = {2024,day_of_year(tp_utc)};
    check_close(eo.distance(),1.47102016e11,0.01_pct);
    check_close(eo.declination()*to_degrees,-23.03,0.1_pct);
 
    // Bergen
    double latitude = 60.3925*to_radians;
    double longitude = 5.323333*to_radians; 
    sun_position sp(tp_utc,latitude ,longitude);
    check_close(sp.zenith_angle()*to_degrees,(90+48.45),1_pct);
    check_close(sp.azimuth_angle()*to_degrees,322.5-180,1_pct);
  } end_test_case()

  begin_test_case(sun_position_test_B) {
    double to_radians = std::numbers::pi/180;
    double to_degrees = 180/std::numbers::pi;
    double earth_tilt = 23.4;
 
    // North pole at northern summer solstice
    time_point tp_utc = {2023,6,21,3,57,0};
    double latitude = 90*to_radians;
    double longitude = 0*to_radians; 
    sun_position sp(tp_utc,latitude ,longitude);
    check_close(sp.zenith_angle()*to_degrees,90-earth_tilt,0.1_pct);

    // South pole at northern winter solstice
    tp_utc = {2023,12,22,3,27,0};
    latitude = -90*to_radians;
    longitude = 0*to_radians; 
    sp = sun_position{tp_utc,latitude ,longitude};
    check_close(sp.zenith_angle()*to_degrees,90-earth_tilt,0.1_pct);
  } end_test_case()
}
