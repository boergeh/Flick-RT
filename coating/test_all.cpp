#include "../environment/unit_test.hpp"
#include "coating_test.hpp"

int main() {
  using namespace flick;
  unit_test t("coating");
  t.include<coating_test>("coating_test");
  t.run_test_cases();
  return 0;
}
