#include "configuration.hpp"

namespace flick {
  begin_test_case(configuration_test_A) {
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
    c_A1.add<std::string>("An_extra_variable","Should not add problem for reading");
    check(c_A1.exists("age"));
    check(not c_A1.exists("not_there"));
    check_throw(c_A1.set<double>("not_there",1));
    my_config_A c_A2;
    
    check_close(c_A1.get<double>("age"), c_A2.get<double>("age"));
    bool stream = true;
    if (stream) {
      write(c_A1,"./tmp.txt");
      c_A2 = read<my_config_A>("./tmp.txt");
      check_close(c_A2.get<double>("age"),49);
      check_close(c_A2.get<double>("shoe"),42);
    }    
  } end_test_case()

  begin_test_case(configuration_test_B) {
    basic_configuration c;
    c.add<double>("list",{0,1,2},"list of numbers");
    check_close(c.get_vector<double>("list").at(1),1);
    std::stringstream ss("# list of numbers ## list = 3 4");
    ss >> c;
    check_close(c.get_vector<double>("list").at(1),4);
  } end_test_case()
}
