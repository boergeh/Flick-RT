#include "../environment/unit_test.hpp"
#include "input.hpp"

namespace flick {
  begin_test_case(input_test) {
    //accurt::api_config c;
    //std::cout <<c;
    //write(accurt::api_config(),"./tmp.txt");
    //auto c = read<accurt::api_config>("./tmp.txt");
    //accurt::input().write("tmp.txt");
    write(accurt::transmittance_api(),"./transmittance_api.txt");
    std::cout << accurt::run(accurt::transmittance_config());
  } end_test_case()
}
