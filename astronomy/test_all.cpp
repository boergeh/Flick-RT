#include "time_point_test.hpp"
#include "sun_position_test.hpp"

int main() {
  using namespace flick;
  unit_test t("astronomy");
  t.include<time_point_test>();
  t.include<sun_position_test_A>();
  t.include<sun_position_test_B>();
  t.run_test_cases();
  return 0;
} 
