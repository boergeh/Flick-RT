#ifndef flick_basic_command
#define flick_basic_command

#include <string>
#include <vector>
#include "../../environment/input_output.hpp"
#include "../../environment/exception.hpp"

namespace flick {
  namespace command {
    class basic_command {
      std::string name_;
      std::vector<std::string> options_;
      std::vector<std::string> arguments_;
    public:
      basic_command(std::string name) : name_{name}{}
      std::string a(size_t n) const {
	if (n >= arguments_.size())
	  return "";
	return arguments_.at(n);
      }
      std::string opt(size_t n) const {
	if (n >= options_.size())
	  return "";
	return options_.at(n);
      }
      void set_arguments(const std::vector<std::string>& input) {
	options_ = get_options(input);
	arguments_ = get_arguments(input);
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
    private:
      using strings = std::vector<std::string>;
      strings get_options(const strings& input) {
	strings options;
	for (size_t i=0; i<input.size(); ++i) {
	  const std::string &o = input.at(i);
	  if (o.substr(0,2) == "--") {
	    options.push_back(o.substr(o.find("=")+1));
	  }
	}
	return options;
      }
      strings get_arguments(const strings& input) {
	strings arguments;
	for (size_t i=0; i<input.size(); ++i) {
	  if (input.at(i).substr(0,2) != "--") {
	    arguments.push_back(input.at(i));
	  }
	}
	return arguments;
      }
    };
  }
  
  template<typename Command>
  bool run(const std::vector<std::string>& args) {
    Command c;
    c.set_arguments(args);
    if (!args.empty() && c.has_name(args.at(0))) {
      c.run();
      std::cout << std::endl;
      return true;
    }
    return false;
  }
}

#endif
