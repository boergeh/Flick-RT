#include "accurt_test.hpp"

int main() {
  using namespace flick;
  unit_test t("accurt"); 
  t.include<accurt_test_A>();
  t.include<accurt_test_B>();
  t.include<accurt_test_C>();
  t.include<accurt_test_D>();
  t.include<accurt_test_E>(); 
  t.include<accurt_test_F>();
  t.include<accurt_test_G>();
  t.run_test_cases();
  return 0;
} 
