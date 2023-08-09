#include "../environment/unit_test.hpp"
#include "parameterized_monodispersed_mie_test.hpp"
#include "monodispersed_mie_test.hpp"
#include "polydispersed_mie_test.hpp"

int main() {
  using namespace flick;
  unit_test t("mie");

  t.include<p_mono_mie_test_A>();
  t.include<p_mono_mie_test_B>();
  
  t.include<mono_mie_bessel_test_A>();
  t.include<mono_mie_bessel_test_B>();
  t.include<mono_mie_test_A>();
  t.include<mono_mie_test_B>();
  t.include<mono_mie_test_C>();
  t.include<mono_mie_test_D>();
  t.include<mono_mie_test_E>();
  t.include<mono_mie_test_F>();
  t.include<mono_mie_test_G>();

  t.include<poly_mie_test_A>();
  t.include<poly_mie_test_B>();
  t.include<poly_mie_test_C>();
  t.include<poly_mie_test_D>();

  t.include<poly_mie_test_t_matrix>();
  t.include<poly_mie_test_no_absorption>();
  
  t.run_test_cases();
  return 0;
} 
