#include "../environment/unit_test.hpp"
#include "radiator_test.hpp"
#include "filter_test.hpp"

int main() {
  using namespace flick;
  unit_test t("radiator");
  t.include<planck_test>();
  t.include<toa_solar_test>();
  t.include<filter_test_A>();
  t.include<filter_test_B>();
  t.include<filter_test_C>();
  t.run_test_cases();
  return 0;
}
