#include "../environment/unit_test.hpp"
#include "mie_test.hpp"

int main() {
  using namespace flick;
  unit_test t("mie");
  t.include<mie_test_A>("mie_test_A");
  t.include<mie_test_B>("mie_test_B");
  t.run_test_cases();
  return 0;
} 
