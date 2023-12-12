#ifndef flick_search_and_replace
#define flick_search_and_replace

#include "input_output.hpp"

namespace flick {
  class parameter_text {
    std::string s_;
    std::string identifier = "#";
  public:
    parameter_text() = default;
    parameter_text(const std::string& s) : s_{s} {}
    std::string get(const std::string& parameter) {
      size_t n = s_.find(parameter+" = ");
      if (n == std::string::npos or parameter.length()==0)
	throw std::runtime_error("\""+parameter+" = \""+" not found in parameter text");
      size_t n_first = n + parameter.length() + 2; 
      size_t n_last = s_.find(identifier, n_first);
      return trim(s_.substr(n_first, n_last-n_first));
    }
    void set(const std::string& parameter, const std::string& value) {
      //std::string p = get(parameter);
      std::string old_str = parameter + " = " + get(parameter);
      std::string new_str = parameter + " = " + value;
      s_.replace(s_.find(old_str),old_str.length(), new_str);
    }
  private:
    std::string trim(std::string s) const {
      const char* t = " \t\n\r\f\v";
      s.erase(0, s.find_first_not_of(t));
      s.erase(s.find_last_not_of(t) + 1);
      return s;
    }
    friend std::ostream& operator<<(std::ostream &os, const parameter_text& t) {
      os << t.s_;
      return os;
    }
    friend std::istream& operator>>(std::istream &is, parameter_text& t) {
      t.s_ = std::string(std::istreambuf_iterator<char>(is), {});
      return is;
    }
  };
}

#endif
