#include "../environment/unit_test.hpp"
#include "emitter_test.hpp"
#include "receiver_test.hpp"
#include "content_test.hpp"

int main() {
  using namespace flick;
  unit_test t("component");
  t.include<emitter_test>();
  t.include<receiver_test>();
  t.include<content_test>();
  t.run_test_cases();
  return 0;
}
