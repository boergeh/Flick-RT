#include "../../environment/unit_test.hpp"
#include "aerosols_test.hpp"

int main() {
  using namespace flick;
  unit_test t("aerosols");
  t.include<aerosols_test>("aerosols_test");
  t.run_test_cases();
  return 0;
}
