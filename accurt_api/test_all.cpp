#include "input_test.hpp"

int main() {
  using namespace flick;
  unit_test t("input");
  t.include<input_test>();
  t.run_test_cases();
  return 0;
} 
