#include "../environment/unit_test.hpp"
#include "ordinary_mc_test.hpp"

int main() {
  using namespace flick;
  unit_test t("transporter");
  t.include<ordinary_mc_test>("ordinary_mc_test");
  t.run_test_cases();
  return 0;
}
