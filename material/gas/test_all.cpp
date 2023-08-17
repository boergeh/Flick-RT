#include "../../environment/unit_test.hpp"
#include "lines_test.hpp"
#include "profile_test.hpp"
#include "atmospheric_state_test.hpp"
#include "air_test.hpp"

int main() {
  using namespace flick;
  unit_test t("gas");
  t.include<atmospheric_state_test>();
  t.include<air_test_o2>();
  t.include<air_test_o3>();
  t.include<profile_test>();
  t.include<lines_test>();
  t.run_test_cases();
  return 0;
}
