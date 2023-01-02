#ifndef flick_plane_parallel
#define flick_plane_parallel

//#include "../numeric/named_bounded_types.hpp"
#include "../transporter/ordinary_mc.hpp"
//#include "distribution.hpp"

namespace flick {
namespace model {
  class layer {
    thickness t_{1};
    std::shared_ptr<material::base> material_;
    std::shared_ptr<coating::base> coating_;
    std::string name_;
    bool active_receivers_{false};
  public:
    layer() = default;
    layer(thickness t, const std::string& name="")
      : t_{t}, name_{name} {
    }
    template <class Material, class... Args>
    void fill(Args... a) {
      material_ = std::make_shared<Material>(a...);
    }
    template <class Coating, class... Args>
    void coat(Args... a) {
      coating_ = std::make_shared<Coating>(a...);
    }
    void activate_receivers() {
      active_receivers_ = true;
    }
    bool active_receivers() const {
      return active_receivers_;
    }
    void name(const std::string& name) {
      name_ = name;
    }
    const std::string& name() const {
      return name_;
    }
    double thickness() const {
      return t_();
    }
    std::shared_ptr<material::base> material() const {
      return material_;
    }
    std::shared_ptr<coating::base> coating() const {
      return coating_;
    }
    friend std::ostream& operator<<(std::ostream &os, const layer& l) {
      os << "\n name            : " << l.name_;
      os << "\n thickness       : " << l.t_();
      os << "\n active receivers: " << (l.active_receivers_? "yes": "no");
      if (l.coating_ != nullptr) {
	os << "\n coating         :" << *l.coating_;
      }
      if (l.material_ != nullptr) {
	os << "\n material        :" << *l.material_;
      }
      return os;
    }
  };

  template <class Coating, class... Args>
  layer bottom_layer(Args... a) {
    model::layer l(thickness{0});
    l.name("bottom");
    l.coat<Coating>(a...);
    return l;
  }

  template<class V>
  class layered_structure {
    V volume_;
    std::shared_ptr<transporter::ordinary_mc> omc_;
  protected:
    std::vector<layer> layers_;
  public:
    layered_structure() = default;
    layered_structure(const layer& bottom)
    {
      V v;
      add_layer(bottom, v);
    }
    const V& volume() const {
      return volume_;
    }
    void transport_radiation(const emitter& em,
			     const std::string& volume_name) {
      if (volume_name == "bottom")
	throw std::runtime_error("Cannot emit from bottom");
      omc_ = std::make_shared<transporter::ordinary_mc>(volume_);
      double g = 0.7;
      omc_->transport_radiation(em,volume_name,g);
    }
    receiver& outward_receiver(const std::string& layer_name) {
      return omc_->outward_receiver(layer_name); 
    }
    receiver& inward_receiver(const std::string& layer_name) {
      return omc_->inward_receiver(layer_name); 
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const layered_structure<V>& s) {
      os << "\n";
      for (size_t i=0; i<s.layers_.size(); ++i) {
	os << "--------------------------------------------------------------";
	os << s.layers_[s.layers_.size()-i-1] << "\n";
      }
      os << "\n";
      os << "Volume three structure:";
      os << "\n" << s.volume();
      return os;
    }
  protected:
    void add_layer(const layer& l, V& v) {
      v.name(l.name());
      if (l.active_receivers()) {
	v.content().inward_receiver().activate();
	v.content().outward_receiver().activate();
      }
      v.content().coat(l.coating());
      v.content().fill(l.material());
      if (layers_.size()>0)
	v.insert(volume_);
      volume_ = v;
      layers_.push_back(l);
    }
  };
   
  class plane_parallel_structure : public layered_structure<semi_infinite_box> {
  public:
    using layered_structure::layered_structure;
    void add_on_top(const layer& l) {
      semi_infinite_box b;  
      pose p = volume().placement();
      p.move_by({0,0,l.thickness()});
      b.move_by(p.position());
      add_layer(l,b);
    }
  };

  //class spherical_structure : public structure<sphere> {
  //  double radius_;
  //public:
  //  void add_layer_outside() {
  //  }
  //};
  //class cylindrical_geometry : structure<cylinder> {
  //  double radius_;
  //public:
  //  void add_layer_on_top() {
  //  }
  // };
}
}

#endif
