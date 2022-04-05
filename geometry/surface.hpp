#ifndef flick_surface
#define flick_surface

#include <iostream>
#include "../numeric/vector.hpp"
#include "../numeric/pose.hpp"

namespace flick {
namespace geometry {
  namespace surface {
    class base {
    public:
      base(){}
      virtual ~base(){}
      pose intersection()
      // The returned pose is the closest intersection point on
      // observer's z_axis (in global coordinates) and the quaternion
      // needed to rotate the global z-axis into the direction of the
      // surface normal at the intersection point. The returned pose
      // thus describes a local coordinate system at the intersection
      // point.
      {
	return intersection_;
      }
      double distance_to_intersection() const {
	return norm(intersection_.position()-observer_.position());
      }
      bool has_intersection() const {
	return has_intersection_;
      }
      bool encloses_observer() const {
	return encloses_observer_;
      }
      virtual void set_observer(const pose& o) = 0;
    protected:
      pose observer_;
      bool has_intersection_{false};
      bool encloses_observer_{false};
      pose intersection_;
    };

    class plane : public base
    // xy-plane at z=0
    {
    public:
      void set_observer(const pose& o) {
	observer_ = o;
	unit_vector normal = {0,0,1};
	double dp = dot(normal,o.z_direction());
	if ((dp > 0 && o.position().z() < 0) || (dp < 0 && o.position().z() > 0)) {
	  has_intersection_ = true;
	  intersection_.move_to(o.position()-o.z_direction()*o.position().z()/dp);
	} else
	  has_intersection_ = false;
	if (o.position().z() < 0)
	  encloses_observer_ = true;
	else
	  encloses_observer_ = false;
      }
    };
    
    class sphere : public base
    // Sphere with radius r centered in origin.
    // See Wikipedia line-sphere instersection.
    {
      double r_{1};
    public:
      sphere(double r) : r_{r} {}
      //static std::shared_ptr<surface::sphere> ptr(double r) {
      //	 return std::make_shared<surface::sphere>(r);
      // }
      void set_observer(const pose& o) {
	observer_ = o;
	const unit_vector &l = o.z_direction();
	const vector &p = o.position();
	double dotlp = dot(l,p);
	double del = pow(dotlp,2) - dot(p,p) + pow(r_,2);

	if (p.r() < r_)
	  encloses_observer_ = true;
	else
	  encloses_observer_ = false;

	if (encloses_observer_)
	  has_intersection_ = true;
	else if (del >= 0 && dotlp < 0)
	  has_intersection_ = true;
	else
	  has_intersection_ = false;

	if (has_intersection_) {
	  double d = -dotlp - sqrt(del); // shortest distance
	  intersection_.move_to(p + l*d);
	  if (encloses_observer_ && dot(intersection_.position(),l) < 0) {
	    d = -dotlp + sqrt(del);
	    intersection_.move_to(p + l*d);
	  }
	  intersection_.rotate_to(normalize(intersection_.position()));
	}
      }
    };
  
  }
  template <class Surface, class... Args>
  std::shared_ptr<Surface> make_surface(Args... a) {
    return std::make_shared<Surface>(a...);
  }
}
}
#endif
