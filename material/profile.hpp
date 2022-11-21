#ifndef flick_profile
#define flick_profile

#include "../numeric/constants.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/pose.hpp"

namespace flick {
  namespace profile {
    class base
    // absorption or scatteringing coefficient profiles 
    {
    protected:
      pe_function vertical_profile_;
      //pe_function slant_profile_;
      vector position_;
      unit_vector direction_;
    public:
      base(const pe_function& vertical_profile, const pose& start)
	: vertical_profile_{vertical_profile},
	  position_{start.position()},
	  direction_{start.direction()} {
      }
      virtual double optical_depth(double distance) const = 0;
      virtual double distance(double optical_depth) const = 0;      
    };
    
    class constant : public base {
    public:
      constant(double profile_value, const pose& start)
	: base{pe_function{profile_value}, start} {
      }
      double optical_depth(double distance) const {
	return vertical_profile_.value() * distance;
      }
      double distance(double optical_depth) const {
	double v = vertical_profile_.value();
	if (v > 0)
	  return optical_depth / v;
	return std::numeric_limits<double>::max(); 
      }
    };
   
    class z_varying : public base {
    public:
      using base::base;
      double optical_depth(double distance) const {
	double mu = zk();
	if (fabs(mu) < std::numeric_limits<double>::epsilon()*10) {
	  double v = vertical_profile_.value(position_.z());
	  return profile::constant(v,pose(position_,direction_)).
	    optical_depth(distance);
	};
	double tau = vertical_profile_.integral(next_z(0),next_z(distance))/mu;
	return tau;
      }
      double distance(double optical_depth) const {
	double mu = zk();
	if (fabs(mu) < std::numeric_limits<double>::epsilon()*10) {
	  double v = vertical_profile_.value(position_.z());
	  return profile::constant(v,pose(position_,direction_)).
	    distance(optical_depth);
	};
	std::optional<double> z =
	  vertical_profile_.integral_limit_b(next_z(0),optical_depth * mu);
	if (z.has_value()) {
	  double zr0 = dot(unit_vector{0,0,1},position_);
	  return (*z - zr0) / mu;
	}
	return std::numeric_limits<double>::max(); 
      }
    private:
      vector next_position(double distance) const {
	return position_ + direction_ * distance;
      }
      double next_z(double distance) const {
	return dot(next_position(distance), unit_vector{0,0,1});
      }
      double zk() const {
	return dot(unit_vector{0,0,1}, direction_);
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
}

#endif
