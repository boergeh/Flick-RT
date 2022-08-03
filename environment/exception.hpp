#ifndef flick_exception
#define flick_exception

namespace flick {
  class exception : public std::exception {
    std::string message_;
  public:
    exception(const std::string& m) : message_{m}{};
    const char* what() const throw ()
    {
      return message_.c_str();
    }
  };
}

#endif
  
