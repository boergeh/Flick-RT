#ifndef flick_rotation
#define flick_rotation

#include "vector.hpp"
#include "quaternion.hpp"

namespace flick {
  vector rotate(const vector& v, const quaternion& rotation) {
    auto p = quaternion{0,v.x(),v.y(),v.z()};
    p = rotation*p*inv(rotation);
    return vector{p.b,p.c,p.d};  
  }
  quaternion rotation_about(const unit_vector& uv, double angle) {
    vector v = uv*sin(angle/2);
    return quaternion{cos(angle/2), v.x(), v.y(), v.z()};
  }      
  quaternion no_rotation() {
    return quaternion{1,0,0,0};
  }
  quaternion rotation_about_x(double theta) {
    return rotation_about({1,0,0},theta);
  }
  quaternion rotation_about_y(double theta) {
    return rotation_about({0,1,0},theta);
  }
  quaternion rotation_about_z(double theta) {
    return rotation_about({0,0,1},theta);
  }
  quaternion rotation_to(const unit_vector& direction) {
    quaternion qz = rotation_about_z(direction.phi());
    unit_vector y_rot = normalize(rotate({0,1,0},qz));
    quaternion qy = rotation_about(y_rot,direction.theta());
    return qy*qz;
  }   

 class rotation {
    quaternion q_{1,0,0,0};
  public:
    rotation() = default;
    rotation(const quaternion& q) : q_{q} {}
    rotation(const unit_vector& direction)
      : q_{rotation_to(direction)} {}
    const quaternion& operator()() const {
      return q_;
    }
    vector apply_rotation_to(const vector& v) const {
      return rotate(v,q_);
    }
    unit_vector x_direction() const {
      return normalize(apply_rotation_to({1,0,0}));
    }
    unit_vector y_direction() const {
      return normalize(apply_rotation_to({0,1,0}));
    }
    unit_vector z_direction() const {
      return normalize(apply_rotation_to({0,0,1}));
    }
    rotation& rotate_to(const unit_vector& uv) {
      q_ = rotation_to(uv);
      return *this;
    }
    rotation& rotate_by(const quaternion& q) {
      q_ = q * q_;
      return *this;
    }
    rotation& rotate_about(const unit_vector& uv, double angle) {
      return rotate_by(rotation_about(uv,angle));
    }
    rotation& rotate_about_local_x(double angle) {
      return rotate_by(rotation_about(x_direction(),angle));
    }
    rotation& rotate_about_local_y(double angle) {
      return rotate_by(rotation_about(y_direction(),angle));
    }
    rotation& rotate_about_local_z(double angle) {
      return rotate_by(rotation_about(z_direction(),angle));
    }
  };
}

#endif
