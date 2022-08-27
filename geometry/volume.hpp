#ifndef flick_volume
#define flick_volume

#include <stdexcept>
#include "boundary.hpp"

namespace flick {
namespace geometry {
  template<class T>
  class volume
  // A volume is confind by its boundary, and can have several other
  // volumes inside, which again can have volumes inside, forming a
  // tree structure
  {
    std::string name_;
    T content_;
    boundary boundary_;
    volume<T>* outer_volume_; 
    std::vector<volume<T>> inner_volumes_;
  public:
    void clear() {
      inner_volumes_.clear();
      inner_volumes_.shrink_to_fit();
    }
    std::string name() const {
      return name_;
    }
    double small_step() {
      return boundary_.small_step();
    }
    pose placement() const { 
      return boundary_.placement();
    }
    size_t n_inner_volumes() const {
      return inner_volumes_.size();
    }
    volume& inner_volume(size_t n) {
      if (n >= inner_volumes_.size())
	throw std::runtime_error("No inner volumes");
      return inner_volumes_[n];
    }
    bool has_outer_volume() const {
      if (outer_volume_ == NULL)
	return false;
      return true;
    }
    volume& outer_volume() const {
      if (!has_outer_volume())
      	throw std::runtime_error("No outer volume");
      return *outer_volume_;
    }
    std::optional<size_t> closest_inner_volume(const pose& observer) const
    // Closest inner volume number (if any) intersecting observer's
    // z-axis.
    {
      std::optional<size_t> n{};
      double d_min = std::numeric_limits<double>::max();
      for (size_t i=0; i < inner_volumes_.size(); ++i) {
	std::optional<pose> p = inner_volumes_[i].boundary_.
	  intersection(observer);
	if (p.has_value()) {
	  double d = norm((*p).position()-observer.position()); 
	  if (d < d_min) {
	    d_min = d;
	    n = i;
	  }
	}
      }
      return n;
    }
    std::optional<pose> intersection(const pose& observer) const
    // Intersection between observer's z-axis and volume
    // boundary. Returned pose includes needed rotation from global
    // z-axis direction to surface normal at intersection point.
    {
      std::optional<size_t> n = closest_inner_volume(observer);
      if (n.has_value()) {
	return inner_volumes_[*n].boundary_.intersection(observer);
      }
      return boundary_.intersection(observer);      
    }
    uniform_intersections get_uniform_intersections(size_t n_reflections) const{
      double cs = boundary_.characteristic_size();
      return uniform_intersections(boundary_, n_reflections,
				   limits{cs*0.1,cs*100});
    }
    volume& name(const std::string name) {
      name_ = name;
      return *this;
    }
    volume& insert(const volume& v) 
    {
      inner_volumes_.emplace_back(v);
      return *this;
    }    
    T& content() {
      return content_;
    }   
    T& operator()() {
      return content_;
    } 
    volume& move_by(const vector& v) {
      boundary_.move_by(v);
      for (size_t i=0; i<inner_volumes_.size(); ++i) {
	inner_volumes_[i].move_by(v);
      }
      return *this;
    }
    volume& rotate_by(const quaternion& rotation,
		      const vector& rotation_center={0,0,0}) {
      boundary_.rotate_by(rotation, rotation_center);
      for (size_t i=0; i<inner_volumes_.size(); ++i) {
	inner_volumes_[i].rotate_by(rotation, rotation_center);
      }
      return *this;
    }
    void set_outer_volume_pointers()
    // To be used with navigator. Note that any outer_volume pointers
    // set inside volume will be lost upon copying a volume object.
    {
      for (size_t i=0; i < inner_volumes_.size(); ++i) {
	inner_volumes_[i].outer_volume_ = this;
	inner_volumes_[i].set_outer_volume_pointers();
      }      
    }
    friend std::ostream& operator<<(std::ostream &os, const volume<T>& v) {
      os << "(0,0) ";
      v.write(os);
      return os;
    }
  protected:
    volume(const boundary& b, const std::string& name) : 
      boundary_{b}, name_{name}, outer_volume_{NULL} {
    }
  private:
    void write(std::ostream& os=std::cout, size_t tree_depth=0) const {
      os << name_ << ", " << inner_volumes_.size() << ", " << boundary_;
      os << '\n';
      for (size_t i=0; i<inner_volumes_.size(); ++i) {
	for (size_t j=0; j<=tree_depth; ++j)
	  os << "  ";
	os << "("<< tree_depth+1 <<","<< i << ") ";
	inner_volumes_[i].write(os,tree_depth+1);
      }
    }
  };

  template<class T>
  class cube : public volume<T> {
  public:
    cube(double side) :
      volume<T>(cubical_boundary(side),"cube") {
    }
  };
  template<class T>
  class sphere : public volume<T> {
  public:
    sphere(double r) :
      volume<T>(spherical_boundary(r),"sphere") {
    }
  };
  template<class T>
  class semi_infinite_box : public volume<T> {
  public:
    semi_infinite_box() :
      volume<T>(boundary{make_surface<surface::plane>()},"semi_infinite_box") {
    }
  };
 
  template<class T>
  class navigator {
    volume<T>* v_;
  public:
    navigator()=default;
    navigator(volume<T> &v) : v_{&v} {
      v_->set_outer_volume_pointers();
    }        
    volume<T>& go_to_outermost_volume() {
      while(v_->has_outer_volume())
	go_outward();
      return *v_;
    }
    volume<T>& go_inward(size_t n=0) {
      v_ = &v_->inner_volume(n);
      return *v_;
    }
    volume<T>& go_outward() {
      v_ = &v_->outer_volume();
      return *v_;
    }    
    volume<T>& go_to(volume<T>& v) {
      v_ = &v;
      return *v_;
    }
    volume<T>& go_to(const std::string& name) {
      if (v_->name()==name)
	return *v_;
      for (size_t i=0; i < v_->n_inner_volumes(); ++i) {
	go_inward(i);
	go_to(name);
      }
      return *v_;
    }
    volume<T>& find(const std::string& name) {
      volume<T>* v0 = v_;
      volume<T>& v_found = go_to(name);
      v_ = v0;
      return v_found;
    }
    volume<T>& current_volume() const {
      return *v_;
    }
    volume<T>& next_volume(const pose& observer) const {
      std::optional<size_t> n = v_->closest_inner_volume(observer);
      if (n.has_value())
	return v_->inner_volume(*n);
      if (v_->has_outer_volume())
	return v_->outer_volume();
      return current_volume();
    }
    std::optional<pose> next_intersection(const pose& observer) const {
      return v_->intersection(observer);
    }
    bool is_moving_inward(const pose& intersection, const pose& observer) const {
      unit_vector surface_normal = intersection.z_direction();
      if (dot(surface_normal,observer.z_direction()) < 0)
	return true;
      return false;
    }
  };
  
}
}

#endif
