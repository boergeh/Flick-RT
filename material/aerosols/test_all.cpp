#include "../../environment/unit_test.hpp"
#include "../ab_functions.hpp"
#include "aerosols_test.hpp"

int main() {
  using namespace flick;
  unit_test t("aerosols");
  t.include<aerosols_test_A>();
  t.include<aerosols_test_B>();
  t.include<aerosols_test_C>();
  t.run_test_cases();
  return 0;
}
