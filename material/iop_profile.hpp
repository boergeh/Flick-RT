#ifndef flick_iop_profile
#define flick_iop_profile

#include "../numeric/constants.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/function.hpp"
#include "../numeric/pose.hpp"

namespace flick {
  class constant_iop {
    double value_;
  public:
    constant_iop(double value) : value_{value} {}
    double optical_depth(double distance) const {
      return value_ * distance;
    }
    double distance(double optical_depth) const {
      double l = optical_depth / value_; 
      if (std::isfinite(l))
	return l;
      return std::numeric_limits<double>::max(); 
    }
    double value(double height) {
      return value_;
    }
  };

  template<class Function>
  class basic_iop_profile {
  protected:
    Function profile_{0};
    double epsilon_ = std::numeric_limits<double>::epsilon()*10;
  public:
    basic_iop_profile() = default;
    basic_iop_profile(const Function& vertical_profile)
      : profile_{vertical_profile} {}
    virtual double optical_depth(const pose& start, double distance) const = 0;
    virtual double distance(const pose& start, double optical_depth) const = 0;
    const std::vector<double>& heights() const {
      return profile_.x();
    }
    size_t size() const {
      return profile_.size();
    }
    void clear() {
      profile_ = Function();
    }
    size_t low_index_near(double height) const {
      return profile_.low_index_near(height);
    }
    double value(double height) const {
      if (height < profile_.x().front() or height > profile_.x().back())
	return 0;
      return profile_.value(height);
    }
    const std::vector<double>& height_grid() const {
      return profile_.x();
    }
    double integral() const {
      return profile_.integral();
    }
    basic_iop_profile<Function>& add(const basic_iop_profile<Function>& p, const std::vector<double>& heights) {
      if (profile_.size() == 0) {
	profile_ = integral_conservative_add(p.profile_, p.profile_, heights);
	profile_.scale_y(0.5);	
      } else {
	profile_ = integral_conservative_add(profile_, p.profile_, heights);
      }
      if (not std::isfinite(profile_.integral()))
	throw std::runtime_error("Possibly less than two points in iop_profile");
      return *this;
    }
  private:
    friend std::ostream& operator<<(std::ostream &os, const basic_iop_profile<Function>& p) {
      os << p.profile_;
      return os;
    }
  };

  template<class Function>
  class iop_z_profile : public basic_iop_profile<Function> {
    using basic_iop_profile<Function>::epsilon_;
    using basic_iop_profile<Function>::profile_;
    unit_vector z_direction_{0,0,1};
  public:
    using basic_iop_profile<Function>::basic_iop_profile;
    double optical_depth(const pose& start, double distance) const {
      double mu = zk(start);
      if (fabs(mu) < epsilon_) {
	double v = profile_.value(start.position().z());
	return constant_iop(v).optical_depth(distance);
      }
      return profile_.integral(next_z(start,0),next_z(start,distance))/mu;
    }
    double distance(const pose& start, double optical_depth) const {
      double mu = zk(start);
      if (fabs(mu) < epsilon_) {
	double v = profile_.value(start.position().z());
	return constant_iop(v).distance(optical_depth);
      }
      std::optional<double> z =
	profile_.integral_limit_b(next_z(start,0),optical_depth * mu);
      if (z.has_value()) {
	double zr0 = dot(z_direction_, start.position());
	return (*z - zr0) / mu;
      }
      return std::numeric_limits<double>::max(); 
    }
  private:
    vector next_position(const pose& start, double distance) const {
      return start.position() + start.direction() * distance;
    }
    double next_z(const pose& start, double distance) const {
      return dot(next_position(start, distance), z_direction_);
    }
    double zk(const pose& start) const {
      return dot(z_direction_, start.direction());
    }
  };
    
    /*
    class atmospheric : public base {
      double planet_radius_;
    public:
      planetary(const pose& start_pose, pe_function vertical_profile, double planet_radius)
	: base{start_pose,vertical_profile}, planet_radius{planet_radius} {
	for (size_t i=0; i < sl.size(); ++i) {
	  slant_profile_.push_back(vertical_profile_.y[i]/r*(dot(r0,uv)+sl[i]));
	}
      }
      double optical_depth(double distance) const {
	return slant_profile(p.direction()).integral(0,distance);
      }
      double distance(double optical_depth) const {
	return slant_profile(p.direction()).integral_limit_b(0, optical_depth);
      }
    private:      
      compute_slant_profile() {
	std::vector<double> tl  = trajectory_length();
	const std::vector<double>& h = vertical_profile_.x();
	const std::vector<double>& coefficient = vertical_profile_.y();
	for (size_t i=0; i < tl.size(); ++i) {
	  double r = planet_radius_ + h[i];
	  double c = coefficient[i] / r * (dot(start_position_,uv)+tl[i]);
	  slant_profile_.insert({tl[i],c});
	}
      }
      std::vector<double> trajectory_length() {
	const std::vector<double>& h = vertical_profile_.x();
	std::vector<double> tl;
	for (size_t i=0; i < h.size(); ++i) {
	  const vector& r0v = start_pose.position();
	  const unit_vector& uv = start_pose.direction();
	  double r = planet_radius_ + h[i];
	  double r0uv = dot(r0,uv);
	  double radicand = pow(r0uv,2)+pow(r,2)-dot(r0,r0),2);
	  if (radicand > 0) {
	    double value;	    
	    if (rd > 0)
	      value = -r0uv + sqrt(radicand);
	    else
	      value = -r0uv - sqrt(radicand);
	    assert(value > 0);
	    tl.push_back(value);
	  }
	}
      assert(tl.size()>0);
      return tl;
      }
    };
    */

}

#endif
