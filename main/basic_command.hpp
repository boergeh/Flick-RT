#ifndef flick_basic_command
#define flick_basic_command

#include <string>
#include <vector>
#include "../environment/input_output.hpp"

namespace flick {
  namespace {
    class basic_command {
      std::vector<std::string> arguments_;
      std::string name_;
    public:
      basic_command(std::string name) : name_{name}{}
      std::string a(size_t n) {
	if (n >= arguments_.size())
	  return "";
	return arguments_.at(n);
      }
      void set_arguments(const std::vector<std::string>& args) {
	arguments_ = args;
      }
      bool has_name(const std::string& name) {
	return (name_ == name);
      }
      void error() {
	std::cout << "\nCannot recognize "+name_+
	  " arguments. Try 'flick help "+name_+"'.\n" << std::endl;  
      }
      void show(const std::string& fname) {
	auto t = read<flick::text>("main/commands/"+fname);
	std::cout << t << std::endl;  
      }
    };
  }
}

#endif
