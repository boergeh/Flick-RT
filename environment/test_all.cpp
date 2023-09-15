#include "unit_test.hpp"
#include "configuration_test.hpp"

int main() {
  using namespace flick;
  unit_test t("environment");
  t.include<configuration_test>();
  t.run_test_cases();
  return 0;
} 
