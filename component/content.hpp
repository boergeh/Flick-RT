#ifndef flick_content
#define flick_content

//#include "../material/material.hpp" // base class
#include "../coating/coating.hpp" // base class
#include "receiver.hpp"

namespace flick {
  class content
  {
    std::shared_ptr<coating::base> coating_;
    //std::shared_ptr<material::base> material_;
    bool has_coating_{false};

    receiver receiver_;
  public:
    template <class Coating, class... Args>
    void set_coating(Args... a) {
      coating_ = std::make_shared<Coating>(a...);
      has_coating_ = true;
    }
    const coating::base& coating() const {
      if (!coating_)
	throw std::runtime_error("content::coating");
      return *coating_;
    }
    bool has_coating() const {
      return has_coating_;
    }
    receiver& receiver() {
      return receiver_;
    }
    
    friend std::ostream& operator<<(std::ostream &os, const content& c) {
      //os << c.has_active_receiver_;
      return os;
    }
  };
}
#endif
