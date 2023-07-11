#include "../environment/unit_test.hpp"
#include "coating_test.hpp"

int main() {
  using namespace flick;
  unit_test t("coating");
  t.include<coating_test_A>();
  t.include<coating_test_B>();
  t.run_test_cases();
  return 0;
}
