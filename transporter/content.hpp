#ifndef flick_content
#define flick_content

//#include "../material/material.hpp" // base class
//#include "../coating/coating.hpp" // base class
#include "emitter.hpp"
#include "receiver.hpp"

namespace flick {  
  class content
  {
    //coating* boundary_coating_;
    //material* filling_material_;
    std::vector<emitter> emitters_;
    std::vector<receiver> receivers_;
    // material_profile mf_;
  public:
    //content& insert(const emitter& em, const vector& relative_position) {
    content& insert(const emitter& em) {
      //em.move_to(relative_position);
      emitters_.emplace_back(em);
      return *this;
    }
    //content& insert(const receiver& re, const vector& relative_position) {
    content& insert(const receiver& re) {
      //re.move_to(relative_position);
      receivers_.emplace_back(re);
      return *this;
    }
    //    content& fill(const material& m, const profile& profile) {
    //   return *this;
    // }
    //content& coat(const coating& c) {
    //  return *this;
    // }
    
    emitter& get_emitter(){return emitters_.at(0);}
    
    
    friend std::ostream& operator<<(std::ostream &os, const content& c) {
      os << c.emitters_.size() << " " << c.receivers_.size();
      return os;
    }
  };
}
#endif
