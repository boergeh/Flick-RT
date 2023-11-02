#include "../environment/unit_test.hpp"
#include "direction_generator_test.hpp"
#include "vector_test.hpp"
#include "histogram_test.hpp"
#include "pose_test.hpp"
#include "range_test.hpp"
#include "sorted_vector_test.hpp"
#include "function_test.hpp"
#include "physics_function_test.hpp"
#include "table_test.hpp"
#include "flist_test.hpp"
#include "distribution_test.hpp"
#include "value_collection_test.hpp"
//#include "mks_units_test.hpp"

int main() {
  using namespace flick;
  unit_test t("numeric");
  //t.include<mks_units_test>("mks_units_test");
  t.include<sorted_vector_test>();
  t.include<function_test_A>();
  t.include<function_test_B>();
  t.include<function_test_C>();
  t.include<function_test_D>();
  t.include<function_test_E>();
  t.include<direction_generator_test>();
  t.include<vector_test>();
  t.include<histogram_test>();
  t.include<pose_test>();
  t.include<range_test>();
  t.include<physics_function_test>();
  t.include<table_test>();
  t.include<flist_test>();
  t.include<distribution_test>();
  t.include<value_collection_test>();
  t.run_test_cases();
  return 0;
} 
