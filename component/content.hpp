#ifndef flick_content
#define flick_content

//#include "../material/material.hpp" // base class
//#include "../coating/coating.hpp" // base class
//#include "emitter.hpp"
#include "receiver.hpp"

namespace flick {
  class content
  {
    //coating::base* boundary_coating_;
    //material::base* filling_material_;
    receiver receiver_;
    bool has_active_receiver_{false};
  public:
    //content& insert(const emitter& em, const vector& relative_position) {
    //void insert(const emitter& em) {
      //em.move_to(relative_position);
    //  emitters_.emplace_back(em);
    //}
    //content& insert(const receiver& re, const vector& relative_position) {
    void activate_receiver() {
      has_active_receiver_ = true; 
    }
    //    content& fill(const material& m, const profile& profile) {
    //   return *this;
    // }
    //content& coat(const coating& c) {
    //  return *this;
    // }
    
    //emitter& get_emitter(){return emitters_.at(0);}
    
    
    friend std::ostream& operator<<(std::ostream &os, const content& c) {
      os << c.has_active_receiver_;
      return os;
    }
  };
}
#endif
