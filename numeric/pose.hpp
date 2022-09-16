#ifndef flick_pose
#define flick_pose

#include "rotation.hpp"

namespace flick {
  class position {
    vector v_{0,0,0};
  public:
    position(const vector& v) : v_{v} {};
    const vector& operator()() const {
      return v_;
    }
    position& move_to(const vector& v) {
      v_ = v;
      return *this;
    }
    position& move_by(const vector& v) {
      v_ += v;
      return *this;
    }
  };
  
  class pose2 {
    position p_{vector{0,0,0}};
    rotation r_{quaternion{1,0,0,0}};
  public:
    pose2() = default;
    pose2(const position& p, const rotation& r) : p_{p},r_{r}{}   
    const position& get_position() const {
      return p_;
    }
    const rotation& get_rotation() const {
      return r_;
    }
    vector apply_rotation_to(const vector& v) const {
      return r_.apply_rotation_to(v);
    }
    pose2 as_observed_by(const pose2& observer) const
    // Do opposite translation and rotation to observer as those
    // applied to reach this pose.
    {
      position p = p_;
      rotation r = r_;
      position p_o = observer.get_position();
      rotation r_o = observer.get_rotation();
      vector v = rotate(p_()-p_o(),inv(r_o()));
      p.move_to(v);
      r.rotate_by(inv(r_o()));
      return pose2(p,r);
    }
    unit_vector x_direction() const {
      return r_.x_direction();
    }
    unit_vector y_direction() const {
      return r_.y_direction();
    }
    unit_vector z_direction() const {
      return r_.z_direction();
    }
    pose2& move_by(const vector& v) {
      p_.move_by(v);
      return *this;
    }
    pose2& move_to(const vector& v) {
      p_.move_to(v);
      return *this;
    }
    pose2& rotate_to(const unit_vector& uv) {
      r_.rotate_to(uv);
      return *this;
    }
    pose2& rotate_by(const quaternion& q) {
      r_.rotate_by(q);
      return *this;
    }
    pose2& rotate_about(const unit_vector& uv, double angle) {
      r_.rotate_about(uv,angle);
      return *this;
    }
    pose2& rotate_about_local_x(double angle) {
      r_.rotate_about_local_x(angle);
      return *this;
    }
    pose2& rotate_about_local_y(double angle) {
      r_.rotate_about_local_y(angle);
      return *this;
    }
    pose2& rotate_about_local_z(double angle) {
      r_.rotate_about_local_z(angle);
      return *this;
    }
  };
  
  class pose {
    // to be replace by pose2
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
