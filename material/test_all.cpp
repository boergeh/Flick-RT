#include "../environment/unit_test.hpp"
#include "material_test.hpp"
#include "iop_profile_test.hpp"
#include "spheres_test.hpp"
#include "normalized_scattering_matrix_fit_test.hpp"
#include "ab_functions_test.hpp"
#include "layered_iops_test.hpp"
#include "z_profile_test.hpp"
#include "mixture_test.hpp"
//#include "earth_atmosphere_test.hpp"

int main() {
  using namespace flick;
  unit_test t("material");
  t.include<material_test_A>();
  t.include<iop_profile_test>();
  t.include<spheres_test_A>();
  t.include<spheres_test_B>();
  t.include<normalized_scattering_matrix_fit_test>();
  t.include<ab_functions_test_A>();
  t.include<ab_functions_test_B>();
  t.include<layered_iops_test_A>();
  t.include<layered_iops_test_B>();
  //t.include<aggregate_test_hg>();
  //t.include<aggregate_test_aerosols>();
  //t.include<aggregate_test_aerosols_hg>();
  t.include<z_profile_test>();
  t.include<mixture_test_aerosols>();
  t.include<mixture_test_hg>();
  t.include<mixture_test_aerosols_hg>();
  t.run_test_cases();
  return 0;
}
