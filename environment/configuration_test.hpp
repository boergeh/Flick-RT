#include "configuration.hpp"

namespace flick {
  begin_test_case(configuration_test) {
    
    enum {age, name};
    struct my_config : public configuration {
      my_config() {
	set<double>(age, 49, "person's age in years");
	
	set<std::string>(name, {"B","Hamre"},
R"(person's name and a very very very very very very very very very very
very long description text.)");
      }
    };
    my_config c;
    write(c,"./tmp.txt");

    my_config c2;
    c2 = read<my_config>("./tmp.txt");
    check_close(c.get<double>(age),c2.get<double>(age));
    //std::cout << std::setprecision(6) << c2;
    
  } end_test_case()
}
