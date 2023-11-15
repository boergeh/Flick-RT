#include "unit_test.hpp"
#include "configuration_test.hpp"
#include "search_and_replace_test.hpp"

int main() {
  using namespace flick;
  unit_test t("environment");
  t.include<configuration_test_A>();
  t.include<configuration_test_B>();
  t.include<search_and_replace_test>();
  t.run_test_cases();
  return 0;
} 
