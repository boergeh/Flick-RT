#include "../../environment/unit_test.hpp"
#include "marine_cdom_test.hpp"

int main() {
  using namespace flick;
  unit_test t("marine_cdom");
  t.include<marine_cdom_test>();
 
  t.run_test_cases();
  return 0;
}
