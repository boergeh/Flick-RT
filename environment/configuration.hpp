#ifndef flick_configuration
#define flick_configuration

#include "input_output.hpp"
#include <map>

namespace flick {
  class basic_parameter {
  public:
    virtual void print(std::ostream &os) = 0;
    virtual void read(std::istream &is) = 0;
    virtual const std::string& description() = 0;
  };

  template<class T>
  class parameter : public basic_parameter {
    std::vector<T> p_;
    std::string description_;
  public:
    parameter(const std::vector<T>& p, const std::string& description="")
      : p_{p}, description_{description} {
    }   
    T p(size_t n) const {
      return p_.at(n);
    }
    size_t size() const {
      return p_.size();
    }
    const std::string& description() {
      return description_;
    }
    void print(std::ostream &os) {
      os << "/* " << description_ << " */\n";
      for (size_t i = 0; i<p_.size(); ++i)
	os << p_[i] << " ";
    }
    void read(std::istream &is) {
      bool is_description = true;
      std::string str;
      while (is_description) {
	is >> str;
	if(str.find("*/") != std::string::npos) {
	  is_description = false;
	}
	if (is.eof())
	  throw std::runtime_error("parameter: missing text qualifier");
      }     
      for (size_t i = 0; i<p_.size(); ++i)
	is >> p_[i];   
    }
  };
  
  class basic_configuration {    
    std::unordered_map<std::string, std::shared_ptr<basic_parameter>> parameters_;
  public:
    template<class T>
    void add(const std::string& name, const std::vector<T>& p, std::string description="") {
      if (empty(description))
	description = parameters_.at(name)->description();
      parameters_[name] = std::make_shared<parameter<T>>(p,description);
    }
    template<class T>
    void add(const std::string& name, const T& p, const std::string& description="") {
      add(name, std::vector<T>{p}, description);  
    }
    template<class T>
    T get(const std::string& name, size_t element=0) const {
      parameter<T> *p = dynamic_cast<parameter<T>*>(&*parameters_.at(name));
      return p->p(element);
    }
    template<class T>
    size_t size(const std::string& name) const {
      parameter<T> *p = dynamic_cast<parameter<T>*>(&*parameters_.at(name));
      return p->size();
    }
    void add_configuration(const basic_configuration& c) {
      for (auto& [name, val] : c.parameters_) {
	std::cout << "adding..";
	parameters_[name] = c.parameters_.at(name);
      }
    }
  private:
    friend std::ostream& operator<<(std::ostream &os,
				    const basic_configuration& c) {
      for (auto& [name, val] : c.parameters_) {
	c.parameters_.at(name)->print(os);
	os << "\n\n";
      }
      return os;
    }   
    friend std::istream& operator>>(std::istream &is,
				    const basic_configuration& c) {
      for (auto& [name, val] : c.parameters_) {
	c.parameters_.at(name)->read(is);
      }
      return is;
    }
  };
}

#endif
