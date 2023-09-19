#include "../environment/unit_test.hpp"
#include "input.hpp"
#include "../material/atmosphere.hpp"

namespace flick {
  begin_test_case(input_test) {
    //accurt::api_config c;
    //std::cout <<c;
    //write(accurt::api_config(),"./tmp.txt");
    //auto c = read<accurt::api_config>("./tmp.txt");
    //accurt::input().write("tmp.txt");
    //write(accurt::transmittance_api(),"./transmittance_api.txt");
    //c_.add_configuration(read<accurt_configuration>("./"+file_name));
    accurt_configuration c;
    auto m = std::make_shared<material::atmosphere>(material::atmosphere::config());
    std::cout << accurt(c,m).relative_irradiance();
  } end_test_case()
}
