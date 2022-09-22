#ifndef flick_content
#define flick_content

#include "../material/material.hpp" 
#include "../coating/coating.hpp" 
#include "receiver.hpp"

namespace flick {
  class content
  {
    std::shared_ptr<coating::base> coating_;
    std::shared_ptr<material::base> material_;
    bool has_coating_{false};
    bool has_material_{false};
    receiver inward_receiver_;
    receiver outward_receiver_;
  public:
    template <class Coating, class... Args>
    void coat(Args... a) {
      coating_ = std::make_shared<Coating>(a...);
      has_coating_ = true;
    }    
    template <class Material, class... Args>
    void fill(Args... a) {
      material_ = std::make_shared<Material>(a...);
      has_material_ = true;
    }
    void fill(std::shared_ptr<material::base> m) {
      material_ = m;
      has_material_ = true;
    }
    coating::base& coating() const {
      return *coating_;
    }
    material::base& material() const {
      return *material_;
    }
    bool has_coating() const {
      return has_coating_;
    }
    bool has_material() const {
      return has_material_;
    }
    receiver& inward_receiver() {
      return inward_receiver_;
    }
    receiver& outward_receiver() {
      return outward_receiver_;
    } 
  };
}

#endif
