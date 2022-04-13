#include "../environment/unit_test.hpp"
#include "planck_test.hpp"

int main() {
  using namespace flick;
  unit_test t("radiator");
  t.include<planck_test>("planck_test");
  t.run_test_cases();
  return 0;
}
