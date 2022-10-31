#include "../environment/unit_test.hpp"
#include "radiator_test.hpp"
#include "filter_test.hpp"

int main() {
  using namespace flick;
  unit_test t("radiator");
  t.include<planck_test>("planck_test");
  t.include<toa_solar_test>("toa_solar_test");
  t.include<filter_test>("filter_test");
  t.run_test_cases();
  return 0;
}
