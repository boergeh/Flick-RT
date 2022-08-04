#include "../environment/unit_test.hpp"
#include "mc_basic_test.hpp"

int main() {
  using namespace flick;
  unit_test t("mc_basic");
  t.include<mc_basic_test>("mc_basic_test");
  t.run_test_cases();
  return 0;
}
