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
    virtual void read(std::istream &is, const std::string& name) = 0;
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
    void read(std::istream &is, const std::string& name) {
      locate_parameter(is, name);
      for (size_t i = 0; i<p_.size(); ++i) {
	is >> p_[i];
      }
    }
  private:
    void skip_description(std::istream &is) const {
      std::string str;
      bool is_description = true;
      while (is_description) {
	is >> str;
	if(str.find(end_qualifier_) != std::string::npos) {
	  is_description = false;
	}
	if (is.eof() or str.empty())
	  throw std::runtime_error("parameter: missing text qualifier");
      }
    }
    void locate_parameter(std::istream &is, const std::string& name) const {
      std::string current_name, equal_sign;
      bool found = false;
      while (not found) {
	skip_description(is);
	is >> current_name;
	if (current_name == name) {
	  found = true;
	}
	is >> equal_sign;
	if (is.eof() or current_name.empty()) {
	  throw std::runtime_error("parameter name '"+name+"' not found");
	}
      }
    }
  };
  
  class basic_configuration {    
    std::map<std::string, std::shared_ptr<basic_parameter>> parameters_;
    std::string begin_qualifier_ = "/*";
    std::string end_qualifier_ = "*/";
    bool unordered_stream_ = false;
  public:
    template<class T>
    void add(const std::string& name, const std::vector<T>& p, std::string description="") {
      if (description.empty() && parameters_.find(name)!=parameters_.end())
	description = parameters_.at(name)->description();
      parameters_[name] = std::make_shared<parameter<T>>(name,p,description);
      parameters_.at(name)->set_text_qualifiers(begin_qualifier_,end_qualifier_);
    }
    template<class T>
    void add(const std::string& name, const T& p, const std::string& description="") {
      add(name, std::vector<T>{p}, description);  
    }
    template<class T>
    void set(const std::string& name, const std::vector<T>& p, const std::string& description="") {
      ensure_exists(name);
      add(name,p,description);
    }
    void set_unordered_stream(bool b) {
      unordered_stream_ = b;
    }
    template<class T>
    void set(const std::string& name, const T& p, const std::string& description="") {
      ensure_exists(name);
      add(name, p, description);  
    }
    template<class T>
    std::vector<T> get_vector(const std::string& name) const {
      ensure_exists(name);
      parameter<T>* p = dynamic_cast<parameter<T>*>(&*parameters_.at(name));
      return p->p();
    }
    bool exists(const std::string& name) const {
      return (parameters_.find(name) != parameters_.end());
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
      begin_qualifier_ = begin;
      end_qualifier_ = end;
    }
  private:
    void ensure_exists(const std::string& name) const {
      ensure(exists(name),"parameter "+name+" not found");
    }
    void ensure(bool b, const std::string& s) const {
      if (not b)
	throw std::runtime_error("configuration " + s);
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const basic_configuration& c) {
      for (auto& [name, val] : c.parameters_) {
	c.parameters_.at(name)->print(os);
	os << "\n\n";
      }
      os << c.begin_qualifier_ <<" end of file "<< c.end_qualifier_;
      return os;
    }   
    friend std::istream& operator>>(std::istream &is,
				    const basic_configuration& c) {
      for (auto& [name, val] : c.parameters_) {
	c.parameters_.at(name)->read(is,name);
	if (c.unordered_stream_) {
	  is.clear();
	  is.seekg(0);
	}
      }
      return is;
    }
  };
}

#endif
