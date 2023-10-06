#ifndef flick_basic_command
#define flick_basic_command

#include <string>
#include <vector>
#include "../../environment/input_output.hpp"
#include "../../environment/exception.hpp"

namespace flick {
  namespace command {
    class basic_command {
      std::vector<std::string> arguments_;
      std::string name_;
    public:
      basic_command(std::string name) : name_{name}{}
      std::string a(size_t n) const {
	if (n >= arguments_.size())
	  return "";
	return arguments_.at(n);
      }
      void set_arguments(const std::vector<std::string>& args) {
	arguments_ = args;
      }
      size_t size() {
	return arguments_.size();
      }
      bool has_name(const std::string& name) const {
	return (name_ == name);
      }
      void error() const {
	std::cout << "\nCannot recognize "+name_+
	  " arguments. Try 'flick help "+name_+"'.\n" << std::endl;  
      }
      void show(const std::string& fname) const {
	auto t = read<flick::text>("main/commands/"+fname);
	std::cout << t << std::endl;  
      }
    };
  }
  
  template<typename Command>
  bool run(const std::vector<std::string>& args) {
    Command c;
    c.set_arguments(args);
    if (!args.empty() && c.has_name(args.at(0))) {
      try {
	c.run();
      } catch (const flick::exception& e) {
	std::cout << "\n Flick exception in " << e.what() << "\n\n";
      } catch (const std::exception& e) {
	std::cout << "\n c++ standard exception: " << e.what() << "\n\n";
      } catch (...) {
	std::cout << "\n Unknown exception" << "\n\n";
      }
      return true;
    }
    return false;
  }
}

#endif
