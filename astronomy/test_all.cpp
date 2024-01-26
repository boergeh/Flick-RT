#include "time_point_test.hpp"
#include "sun_position_test.hpp"

int main() {
  using namespace flick;
  unit_test t("astronomy");
  t.include<time_point_test_A>();
  t.include<time_point_test_B>();
  t.include<time_point_test_C>();
  t.include<time_point_test_D>();
  t.include<time_point_test_E>();
  t.include<sun_position_test_A>();
  t.include<sun_position_test_B>();
  t.run_test_cases();
  return 0;
} 
