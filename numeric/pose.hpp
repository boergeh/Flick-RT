#ifndef flick_pose
#define flick_pose

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
    
  class pose {
    vector position_{0,0,0};
    quaternion rotation_{1,0,0,0};
  public:
    pose() = default;
    pose(const vector& position, const quaternion& rotation)
      : position_{position}, rotation_{rotation} {}
    pose(const vector& position, const unit_vector& direction)
      : position_{position}, rotation_{rotation_to(direction)} {}
    const vector& position() const {return position_;}
    const quaternion& rotation() const {return rotation_;}
    vector apply_pose_rotation_to(const vector& v) const {
      return rotate(v,rotation_);
    }
    pose as_observed_by(const pose& observer) const
    // Do opposite translation and rotation to observer as those
    // applied to reach this pose.
    {
      pose p = *this;
      vector v = rotate(position_- observer.position(),
			inv(observer.rotation()));
      p.move_to(v);
      p.rotate_by(inv(observer.rotation()));
      return p;
    }
    unit_vector x_direction() const {
      return normalize(apply_pose_rotation_to({1,0,0}));
    }
    unit_vector y_direction() const {
      return normalize(apply_pose_rotation_to({0,1,0}));
    }
    unit_vector z_direction() const {
      return normalize(apply_pose_rotation_to({0,0,1}));
    }
    unit_vector direction() const {
      return z_direction();
    }
    pose& move_by(const vector& v) {
      position_ += v;
      return *this;
    }
    pose& move_to(const vector& position) {
      position_ = position;
      return *this;
    }
    pose& rotate_to(const unit_vector& uv) {
      rotation_ = rotation_to(uv);
      return *this;
    }
    pose& rotate_by(const quaternion& rotation) {
      rotation_ = rotation * rotation_;
      return *this;
    }
    pose& rotate_about(const unit_vector& uv, double angle) {
      return rotate_by(rotation_about(uv,angle));
    }
    pose& rotate_about_local_x(double angle) {
      return rotate_by(rotation_about(x_direction(),angle));
    }
    pose& rotate_about_local_y(double angle) {
      return rotate_by(rotation_about(y_direction(),angle));
    }
    pose& rotate_about_local_z(double angle) {
      return rotate_by(rotation_about(z_direction(),angle));
    }
  };
  
  std::ostream& operator<<(std::ostream& ostr, const pose& p) {
    ostr << p.position() << " "<< p.rotation();
    return ostr;
  }
  void show_position_and_directions(const pose& p) {
    std::cout << std::endl;
    std::cout <<  p.position() << std::endl;
    std::cout <<  p.x_direction() << std::endl;
    std::cout <<  p.y_direction() << std::endl;
    std::cout <<  p.z_direction();
  }
  pose global_observer() {return pose();}
  const pose default_pose = pose{};
}

#endif
