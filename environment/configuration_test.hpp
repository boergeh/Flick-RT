#include "configuration.hpp"

namespace flick {
  begin_test_case(configuration_test) {
    struct my_config_A : public basic_configuration {
      my_config_A() {
	add<double>("age", 49, "person's age in years");	
	add<double>("shoe", 42, "person's shoe size");	
	add<std::string>("name", {"B","Hamre"},
			 R"(person's name and a very very very very very very very very very very
	very long description text.)");
      }
    };
    
    struct my_config_B : public basic_configuration {
      my_config_B() {
	add<std::string>("bike", "candondale", "person's bike");	
      }
    };

    struct my_config_C : public basic_configuration {
      my_config_C() {
	add_configuration(my_config_A());
	add_configuration(my_config_B());
      }
    };
        
    my_config_A c_A1;
    write(c_A1,"./tmp.txt");
    my_config_A c_A2;
    c_A2 = read<my_config_A>("./tmp.txt");
    check_close(c_A1.get<double>("age"), c_A2.get<double>("age"));
    //std::cout << std::setprecision(6) << c_A2;
    //std::cout << c_A2.get<double>("shoe");
    //std::cout << my_config_C();
    
  } end_test_case()
}
