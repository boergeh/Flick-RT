#ifndef flick_configuration
#define flick_configuration

#include "input_output.hpp"
#include <map>

namespace flick {
  class basic_parameter {
  protected:
    std::string begin_qualifier_ = "/*";
    std::string end_qualifier_ = "*/";
  public:
    virtual void print(std::ostream &os)=0;
    virtual void read(std::istream &is) = 0;
    virtual const std::string& description() = 0;
    void set_text_qualifiers(const std::string& begin, const std::string& end) {
      begin_qualifier_ = begin;
      end_qualifier_ = end;
    }

  };

  template<class T>
  class parameter : public basic_parameter {
    std::string name_;
    std::vector<T> p_;
    std::string description_;
  public:
    parameter(const std::string& name, const std::vector<T>& p, const std::string& description="")
      : name_{name}, p_{p}, description_{description} {
    }   
    T p(size_t n) const {
      return p_.at(n);
    }
    std::vector<T> p() const {
      return p_;
    }
    size_t size() const {
      return p_.size();
    }
    const std::string& description() {
      return description_;
    }
    void print(std::ostream &os) {
      os << begin_qualifier_ << " " << description_ << " " << end_qualifier_ << "\n";
      os << name_ << " = ";
      for (size_t i = 0; i<p_.size(); ++i)
	os << p_[i] << " ";
    }
    void read(std::istream &is) {
      bool is_description = true;
      std::string str;
      while (is_description) {
	is >> str;
	if(str.find(end_qualifier_) != std::string::npos) {
	  is_description = false;
	}
	if (is.eof())
	  throw std::runtime_error("parameter: missing text qualifier");
      }
      is >> name_ >> str;
      for (size_t i = 0; i<p_.size(); ++i)
	is >> p_[i];   
    }
  };
  
  class basic_configuration {    
    std::map<std::string, std::shared_ptr<basic_parameter>> parameters_;
    std::string qualifier_begin_ = "/*";
    std::string qualifier_end_ = "*/";
  public:
    template<class T>
    void add(const std::string& name, const std::vector<T>& p, std::string description="") {
      if (empty(description) && parameters_.find(name)!=parameters_.end())
	description = parameters_.at(name)->description();
      parameters_[name] = std::make_shared<parameter<T>>(name,p,description);
      parameters_.at(name)->set_text_qualifiers(qualifier_begin_,qualifier_end_);
    }
    template<class T>
    void add(const std::string& name, const T& p, const std::string& description="") {
      add(name, std::vector<T>{p}, description);  
    }
    template<class T>
    void set(const std::string& name, const std::vector<T>& p, std::string description="") {
      add(name,p,description);
    }
    template<class T>
    void set(const std::string& name, const T& p, const std::string& description="") {
      add(name, p, description);  
    }
    template<class T>
    std::vector<T> get_vector(const std::string& name) const {
      parameter<T>* p;
      try {
	p = dynamic_cast<parameter<T>*>(&*parameters_.at(name));
      } catch (const std::exception& e) {
	throw std::runtime_error("configuration parameter "+name+" not found");
      }
      return p->p();
    }
    template<class T>
    T get(const std::string& name, size_t element=0) const {
      std::vector<T> p = get_vector<T>(name);
      return p.at(element);
    }
    template<class T>
    size_t size(const std::string& name) const {
      parameter<T> *p = dynamic_cast<parameter<T>*>(&*parameters_.at(name));
      return p->size();
    }
    void add_configuration(const basic_configuration& c) {
      for (auto& [name, val] : c.parameters_) {
	parameters_[name] = c.parameters_.at(name);
      }
    }
    void set_text_qualifiers(const std::string& begin, const std::string& end) {
      for (auto& [name, val] : parameters_) {
	parameters_.at(name)->set_text_qualifiers(begin,end);
      }
      qualifier_begin_ = begin;
      qualifier_end_ = end;
    }
  private:
    friend std::ostream& operator<<(std::ostream &os,
				    const basic_configuration& c) {
      for (auto& [name, val] : c.parameters_) {
	c.parameters_.at(name)->print(os);
	os << "\n\n";
      }
      os << "#";
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
