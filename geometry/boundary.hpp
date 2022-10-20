#ifndef flick_boundary
#define flick_boundary
#include "surface.hpp"
#include "../numeric/direction_generator.hpp"
#include <optional>

namespace flick {
namespace geometry {
  const bool inside_out{true};  
  class boundary
  // A boundary is built from one or more surfaces. The boundary has
  // its center in the origin until it has been moved.
  {
    struct element {
      std::shared_ptr<surface::base> surface_ptr;
      pose placement;
      bool inside_out{false};
    };
    std::vector<element> elements_;
    pose placement_{{0,0,0},no_rotation()};
    double characteristic_size_{1};
  public:
    boundary() = default;   
    boundary(std::shared_ptr<surface::base> s,
	     const pose& placement=pose{}, bool inside_out=false) {
      add(s,placement,inside_out);
    }
    boundary& set_characteristic_size(double cs) {
      characteristic_size_ = cs;
      return *this;
    }
    double characteristic_size() const {
      return characteristic_size_;
    }
    boundary& add(std::shared_ptr<surface::base> s,
		  const pose& placement=pose{},
		  bool inside_out=false) {
      element e{};
      e.placement = placement;
      e.surface_ptr = s;
      e.inside_out = inside_out;
      elements_.emplace_back(e);
      return *this;
    }
    double small_step() const {
      return characteristic_size()*1e3*std::numeric_limits<double>::epsilon();
    }       
    pose placement() const
    // Get position and orientation of local boundary coordinate
    // system
    {  
      return placement_;
    }
    std::optional<pose> intersection(const pose& observer) const
    // Intersection between observer's z-axis and innermost surface
    // element (if any). Returned pose includes needed rotation from
    // global z-axis direction to surface normal at intersection
    // point.
    {
      std::optional<size_t> surface_n{};
      if (is_enclosed(observer)) {
      	surface_n = closest_surface(observer);
      } else {
	surface_n = enclosed_by_all_others(observer);
      }
      if (surface_n.has_value()) {
	surface_as_observed_by(*surface_n, observer);
	pose global_observer = get_globally_observed_intersection(*surface_n);
	if (elements_.at(*surface_n).inside_out)
	  global_observer.rotate_about_local_x(constants::pi);
	return std::optional<pose>{global_observer};
      }
      return std::optional<pose>{};
    }
    boundary& move_by(const vector& v) {
      for (size_t i=0; i<elements_.size(); ++i)
	elements_[i].placement.move_by(v);
      placement_.move_by(v); 
      return *this;
    }
    boundary& rotate_by(const quaternion& rotation, const vector& rotation_center={0,0,0})
    {
      pose p0 = placement_;
      p0.move_to(rotation_center);
      p0.rotate_by(inv(rotation));
      for (size_t i=0; i<elements_.size(); ++i) {
	elements_[i].placement = elements_[i].placement.as_observed_by(p0);
	elements_[i].placement.move_by(p0.position());
      }
      placement_ = placement_.as_observed_by(p0);
      placement_.move_by(p0.position());
      return *this;
    }
    friend std::ostream& operator<<(std::ostream &os, const boundary& b) {
      os << b.placement_ << ", " << b.elements_.size();
      return os;
    }
  private:
    pose get_globally_observed_intersection(size_t n) const {
      const pose& p0 = elements_.at(n).surface_ptr->intersection();
      const pose& p1 = global_observer().as_observed_by(elements_.at(n).placement);
      return p0.as_observed_by(p1);
    }
    surface::base* surface_as_observed_by(size_t n, const pose& observer) const {
      surface::base *s = elements_.at(n).surface_ptr.get();
      s->set_observer(observer.as_observed_by(elements_.at(n).placement));
      return s;
    }
    std::optional<size_t> closest_surface(const pose& observer) const {
      double distance = std::numeric_limits<double>::max();
      std::optional<size_t> surface_n{};
      for (size_t i=0; i < elements_.size(); ++i) {
	surface::base *s = surface_as_observed_by(i,observer);
	if (s->has_intersection()) {
	  double d = s->distance_to_intersection();
	  if (d < distance) {
	    distance = d;
	    surface_n = i;
	  }
	}
      }
      return surface_n;
    }
    std::optional<size_t> enclosed_by_all_others(const pose &observer) const
    // Returns number of the closest surface if it has intersection
    // point with observer's z-axis and the intersection point is
    // enclosed by all other surfaces in the boundary.
    {
      std::optional<size_t> surface_n{};
      double d_min = std::numeric_limits<double>::max();
      for (size_t i=0; i < elements_.size(); ++i) {
	surface::base *s = surface_as_observed_by(i,observer);
	bool hi = s->has_intersection();
	double d = s->distance_to_intersection();
	if (hi && d < d_min) {
	  pose goi = get_globally_observed_intersection(i);
	  if (is_enclosed(goi,i)) {
	    d_min = d;
	    surface_n = i;
	  }
	}      
      }
      return surface_n;
    }
    bool is_enclosed(const pose& observer,
		     size_t skip_n=std::numeric_limits<size_t>::max()) const {
      for (size_t i=0; i < elements_.size(); ++i) {
	if (i!=skip_n) {
	  surface::base *s = surface_as_observed_by(i,observer);
	  if (!s->encloses_observer() && !elements_.at(i).inside_out) {
	    return false;
	  }
	}
      }
      return true;
    }
  };
  
  class cubical_boundary : public boundary {
  public:
    cubical_boundary(double side) {
      set_characteristic_size(side);
      auto s = make_surface<surface::plane>();
      double d = side/2;
      add(s,{{0,0,d},rotation_to({0,0,1})});
      add(s,{{0,0,-d},rotation_to({0,0,-1})});
      add(s,{{0,d,0},rotation_to({0,1,0})});
      add(s,{{0,-d,0},rotation_to({0,-1,0})});
      add(s,{{d,0,0},rotation_to({1,0,0})});
      add(s,{{-d,0,0},rotation_to({-1,0,0})});
    }
  };

  class spherical_boundary : public boundary {
  public:
    spherical_boundary(double radius) {
      set_characteristic_size(radius);
      add(make_surface<surface::sphere>(radius));
    }
  };

  class uniform_intersections
  // Uniformly distributed boundary intersections as result of placing
  // boundary inside a sphere with white lambertian walls
  { 
    boundary boundary_;
    size_t n_reflections_;
    double outer_sphere_area_;
    std::vector<pose> intersections_;
  public:
    uniform_intersections(const boundary& b, size_t n_reflections,
			  limits outer_sphere_limits=limits{0.1,1000})
      : boundary_{b}, n_reflections_{n_reflections}
    {
      double r = outer_sphere_limits.lower();
      intersections_ = find_intersections(r);
      while (intersections_.size()==0 && r <= outer_sphere_limits.upper()) {
	r *= 2;
	intersections_ = find_intersections(r);
      }
      outer_sphere_area_ = 4*constants::pi*pow(r,2);
    }
    std::vector<pose> get() const {
      return intersections_;
    }
    double enclosed_volume() const
    // Using divergence theorem to relate surface to volume.
    {
      double v = 0;
      size_t n =  intersections_.size();
      for (size_t i = 0; i < n; i++) {
	const pose& p = intersections_.at(i);
	v += dot(p.position(),p.z_direction());
      }
      return 0.5/3*v/n_reflections_*outer_sphere_area_;
    }
    double boundary_area() const {
      return 0.5*intersections_.size()*outer_sphere_area_/n_reflections_;
    }
    vector center_of_gravity() const
    // Using divergence theorem to relate surface to center of gravity.
    {
      double cg_x = 0;
      double cg_y = 0;
      double cg_z = 0;
      for (size_t i = 0; i < intersections_.size(); i++) {
	const pose& p = intersections_.at(i);
	vector vx = 0.5*pow(p.position().x(),2)*vector{1,0,0};
	vector vy = 0.5*pow(p.position().y(),2)*vector{0,1,0};
	vector vz = 0.5*pow(p.position().z(),2)*vector{0,0,1};
	unit_vector n = p.z_direction();
	cg_x += dot(vx,n);
	cg_y += dot(vy,n);
	cg_z += dot(vz,n);
      }
      return 0.5*vector{cg_x,cg_y,cg_z}/
	n_reflections_*outer_sphere_area_/enclosed_volume();      
    }    
    friend std::ostream& operator<<(std::ostream &os, const uniform_intersections& ui) {
      for (size_t i = 0; i<ui.intersections_.size(); i++) {
	const pose& p = ui.intersections_.at(i);
	os << p.position() << "   " << p.z_direction() << '\n';
      }
      return os;
    }
  private:
    std::vector<pose> find_intersections(double radius)
    // Note that uncertainties increase when boundaries are moved away
    // from the origin.
    {
      const boundary& b = boundary_;
      std::vector<pose> intersections;
      direction_generator dg;
      auto sphere = std::make_shared<surface::sphere>(radius);
      boundary reflector{sphere};
      reflector.move_by(b.placement().position());
      pose observer = b.placement().rotate_to(dg.isotropic());
      size_t n = 0;
      while (n < n_reflections_) {
	std::optional<pose> pb = b.intersection(observer);
	std::optional<pose> pr = reflector.intersection(observer);
	if (pb.has_value()) {
	  intersections.emplace_back(*pb);
	  observer.move_to((*pb).position());
	  observer = step_away_from_surface(observer,(*pb).z_direction(),
					    b.small_step());
	} else if (pr.has_value()) {
	  observer = b.placement().rotate_to(dg.isotropic());
	  observer.move_to((*pr).position());
	  observer = step_away_from_surface(observer,(*pr).z_direction(),
					    -b.small_step());	  
	  observer.rotate_to(dg.lambertian(-(*pr).z_direction()));	  
	  n++;
	} else {
	  return std::vector<pose>{};
	}
      }
      return intersections;
    }
    pose step_away_from_surface(pose p, const unit_vector& surface_normal,
				 double step_length) {
      int s = +1;
      double d = dot(p.z_direction(),surface_normal);
      if (d < 0)
	s = -1;
      return p.move_by(s*surface_normal*step_length);
    }
  };
}
}
#endif
