#ifndef flick_input_output
#define flick_input_output

#include <fstream>
#include <iostream>
#include <iomanip>

namespace flick {
  std::string path() {
    std::string path;
    char *p = getenv("FLICK_PATH");
    if (p==NULL)
      throw std::invalid_argument("FLICK_PATH not found. See Prerequisites");
    return std::string(p);
  }

  class text {
    std::string t_;
    friend std::istream& operator>>(std::istream &is, text& t) {
      t.t_= std::string(std::istreambuf_iterator<char>(is),
	std::istreambuf_iterator<char>());
      return is;
    }
    friend std::ostream& operator<<(std::ostream &os,text& t) {
      os << t.t_;
      return os;
    }
  };
  
  template<typename T>
  T read(std::string file_name) {
    if (file_name.substr(0,2)!="./")
      file_name = path()+"/"+file_name;     
    std::ifstream ifs(file_name);
    if (!ifs)
      throw std::invalid_argument(file_name+" not found");
    T t;
    ifs >> t;
    ifs.close();
    return t;	
  }  

  template<typename T>
  void write(T& t, const std::string& file, size_t precision=5) {
    std::ofstream ofs(path()+"/"+file);
    if (!ofs)
      throw std::invalid_argument(file+" could not be opened");
    ofs << std::setprecision(precision) << t;
    ofs.close();
  }
}


#endif
  
