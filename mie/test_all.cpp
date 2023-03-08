#include "../environment/unit_test.hpp"
#include "mie_test.hpp"

int main() {
  using namespace flick;
  unit_test t("mie");
  t.include<mie_test>("mie_test");
  t.run_test_cases();
  return 0;
} 
