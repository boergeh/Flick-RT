#include "../../environment/unit_test.hpp"
#include "command_test.hpp"

int main() {
  using namespace flick;
  unit_test t("command");
  t.include<command_test>();
  t.run_test_cases();
  return 0;
}
