#include "../../environment/unit_test.hpp"
#include "lines_test.hpp"
#include "profile_test.hpp"
#include "atmosphere_state_test.hpp"

int main() {
  using namespace flick;
  unit_test t("gas");
  t.include<atmosphere_state_test>("atmosphere_state_test");
  //t.include<profile_test>("profile_test");
  //t.include<lines_test>("lines_test");
  t.run_test_cases();
  return 0;
}
