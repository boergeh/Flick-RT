#ifndef flick_configuration
#define flick_configuration

#include "input_output.hpp"

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
  
  class configuration {    
    std::vector<std::shared_ptr<basic_parameter>> parameters_;
  public:
    template<class T>
    void set(int n, const std::vector<T>& p, std::string description="") {
      ensure_size(n);
      if (empty(description))
	description = parameters_[n]->description();
      parameters_[n] = std::make_shared<parameter<T>>(p,description);    
    }    
    template<class T>
    void set(int n, const T& p, const std::string& description="") {
      set(n, std::vector<T>{p}, description);  
    }    
    template<class T>
    T get(int n, size_t element=0) const {
      parameter<T> *p = dynamic_cast<parameter<T>*>(&*parameters_.at(n));
      return p->p(element);
    }
  private:
    friend std::ostream& operator<<(std::ostream &os,
				    const configuration& c) {
      for (size_t i = 0; i<c.parameters_.size(); ++i) {
	c.parameters_[i]->print(os);
	os << "\n\n";
      }
      return os;
    }
    void ensure_size(size_t n) {
      if (parameters_.size() <= n)
	parameters_.resize(n+1); 
    }
    friend std::istream& operator>>(std::istream &is,
				    const configuration& c) {
      for (size_t i = 0; i<c.parameters_.size(); ++i) {
	c.parameters_[i]->read(is);
      }
      return is;
    }
  };
}

#endif
