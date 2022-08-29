#ifndef flick_material
#define flick_material

#include "../environment/input_output.hpp"
#include "../numeric/constants.hpp"
#include "../numeric/function.hpp"
#include "../numeric/physics_function.hpp"
#include "../numeric/pose.hpp"
#include "../numeric/bounded_type.hpp"
#include "../polarization/rayleigh_mueller.hpp"

namespace flick {
  using wavelength = bounded_type<double, zero, std::exa>;
namespace material {
  class base {
  protected:
    wavelength wavelength_{500e-9};
    pose pose_;
  public:
    void set(const pose& p) {
      pose_ = p;
    }
    void set(const wavelength& wl) {
      wavelength_ = wl;
    }
    double angle(const unit_vector& scattering_direction) {
      return acos(dot(pose_.z_direction(),scattering_direction));
    }
     //void direction(const unit_vector& d) {
    //  direction_ = d;
    //}
    //virtual double absorption_coefficient() = 0; //remove?
    //virtual double scattering_coefficient() = 0; //remove?
    //coefficient coefficient() {}

    virtual double absorption_optical_depth(double distance) {
      // tbi
      return 0;
    }
    virtual double scattering_optical_depth(double distance) {
      return 0;
    }
    virtual double absorption_distance(double absorption_optical_depth) {
      // tbi
      return std::numeric_limits<double>::max();
    }
    virtual double scattering_distance(double scattering_optical_depth) {
      // tbi
      return std::numeric_limits<double>::max();
    }
    virtual mueller mueller_matrix(const unit_vector& scattering_direction) {
      mueller m;
      m.add(0,0,1/(4*constants::pi));
      return m;
    }
    virtual double sampling_asymmetry_factor() {
      return 0.8;
    }
    /*
    unit_vector proposed_scattering_direction(double random_polar,
					      double random_azimuth) {
      double g = proposed_asymmetry_factor_;
      double polar_angle = henyey_greenstein{g}.inverted_accumulated_angle(random_polar);
      double azimuth_angle = random_azimuth*2*constants::pi;
      return {polar_angle, azimuth_angle};
    }
    */
    /*
    double proposed_polar_angle(double random) {
      double g = proposed_asymmetry_factor_;
      return henyey_greenstein{g}.inverted_accumulated_angle(random);
    }
    double proposed_azimuth_angle(double random) {
      double isotropic = random * 2 * constants::pi;
      return isotropic;
    }
    */
  };
  class vacuum : public base {
  public:
    double refractive_index() {
      return 1;
    }
  };
  /*
  void stream_basic_material(std::ostream &os, basic_material& bm) {
    os << "absorption coefficient " << bm.absorption_coefficient()
       << ", scattering coefficient " << bm.scattering_coefficient();
    //return os;
  }
  */
  /*
  class gas : public material {
    std::vector<double> wavelengths_{500e-9};
  public:
    gas(const pe_function& total_pressure, const pe_function& partial_pressure) {
    }
    //virtual iops iops(wl,sdf) {
    //}
  };
    */  

  
  
  class basic_profile {
  public:
    virtual double column_density(double column_length,
				  const pose& p) = 0;
    virtual std::optional<double> column_length(double column_density,
				 const pose& p) = 0;
  };
  
  class constant_profile : public basic_profile {
    double volumetric_density_;
  public:
    constant_profile(double volumetric_density)
      : volumetric_density_{volumetric_density} {
    }
    double column_density(double column_length, const pose& p) {
      return volumetric_density_ * column_length;
    }
    std::optional<double> column_distance(double column_density, const pose& p) {
      return column_density / volumetric_density_;
    }    
  };  
  /*
  class exponential_profile : public profile {
    function<piecewise_exponential_like> volumetric_density;
  public:
    exponential_profile(const point& p_low, point p_high)
      : volumetric_density{{p_low.x(),p_high.x()},{p_low.y(),p_high.y()}} {
    }
    void append(const point& p) {
      volumetric_density.append(p);
    }
    double column_density(double column_length, const pose& p) {
      double u = p.direction().u();
      double h1 = p.position().z();
      double h2 = h1 + column_length * u;
      if (fabs(u) < 1e3*std::numeric_limits<double>::epsilon())
      	return volumetric_density.value(h1) * column_length;
      return volumetric_density.integral(h1, h2) / u;
    }
    std::optional<double> column_length(double column_density, const pose& p) {
      double u = p.direction().u();
      double h1 = p.position().z();
      if (fabs(u) < 1e3*std::numeric_limits<double>::epsilon())
	return column_density / volumetric_density.value(h1);
      return volumetric_density.integral_limit_b(h1, column_density * u);
    }    
  };
  */
 

  
  /*
  class curved_exponential_profile : public exponential_profile {
    // see plot air mass.m in blue_sky
    vector r_v;
    unit_vector d_uv;
    double r{0};
    vector& r_uv;
    size_t nc{0};
    double mu{0};
    double k{0};
    double d_next{0};
    double m_next{0};    
    double characteristic_length;
    double epsilon;
    
  public:
    curved_exponential_profile(const point& low, const point& high) 
      : characteristic_length{high.height-low.height},
	epsilon{characteristic_length*1e-9}
    {
      add_point(low.height, low.value);
      add_point(low,value_b);
    }      
    double mass(double distance) {
      double m = 0;
      double d = 0;      
      while (true) {
	double d_next = distance_to_next_layer();
	if (d+d_next > distance) {
	  break;
	}
	d += d_next;
	m += mass_in_current_layer(d_next);
	move_to_next_layer(d_next);
      }
      m += mass_in_current_layer(distance-d);
      return m;
    }
    double distance(double mass) {      
      double d = 0;
      double m = 0;
      while (true) {
	d_next = distance_to_next_layer();
	m_next = mass_in_current_layer(d_next);
	if (m + m_next > mass)
	  break;
	m += m_next;
	d += d_next;
	move_to_next_layer(d_next);
      }
      d += distance_in_current_layer(mass-m);
      return d;
    }
    void set_positon_and_direction(const vector& position, const unit_vector& direction) {
      r_v = positon;
      d_uv = direction;
      r = norm(r_v);
      r_uv = normalize(r);
      nc = find_point_below(r);
      mu = dot(r_uv, d_uv);
      k = slope(nc);
    } 
  private:
    size_t find_point_below(double h) {
      size_t n = n_start;
      bool below_bottom = r < points[0].height;
      bool above_top = r > points.back().height;
      if (below_bottom)
	return 0;
      if (above_top)
	return heig.size()-1;
      if (inside_layer(n))
	return n;
      if (in_layer_below(n))
	return n-1;
      if (in_layer_above(n))
	return n+1;
      while (!in_layer(n)) {
	if (r < height[n])
	  n /= 2;
	else
	  n += (heights.size()-n)/2;
      }
      return n;
    }
    bool in_layer_below(size_t n) {
      return (r < points[n].height && r > points[n-1].height);
    }
    bool in_layer_above(size_t n) {
      return (r > points[n+1].height && r < points[n+2].height);
    }
    bool in_layer(size_t n) {
      return (points[n].height < r && points[n+1].height > r);
    }
    double distance_to_next_layer() {
      double dl = distance_to_layer(heights[nc]);
      double dh = distance_to_layer(heights[nc+1]);
      if (dl < dh)
	return dl;
      return dh;
    }
     double distance_to_layer(double layer_radius) {
      double dicriminant = mu-pow(layer_radius/r,2)-1;
      if (discriminant < 0)
	return std::numeric_limits<double>::max();
      double d = r*(sqrt(discriminant)-mu); 
      if (d <= 0)
	return std::numeric_limits<double>::max();
      return d;
    }
    double mass_in_current_layer(double d) {
      double a = k*mu;
      if (fabs(a) < epsilon)
	return density[nc]*d;
      return density[nc]/a*(exp(a*d)-1);
    }    
    void move_to_next_layer(duble d_next) {
      vector r_v_new = r_v + d_uv*(d_next+epsilon)
      set_current_position_and_direction(r_v_new, d_uv);
    }
    double distance_in_current_layer(double m) {
      double a = k*mu;
      if (fabs(a) < epsilon)
	return m/density[m];
      return log(m/denisty[n]*k*mu+1)/a;      
    }        
    double slope(size_t n) {
      if (n >= height_.size())
	n = height_.size()-1;
      return log(density_[n+1]/density_[n])/(height_[n+1]-height_[n]);
    }
  };

  class material : public fabric {
    std::shared_ptr<profile> profile_{std::make_shared<constant_profile>()};
  public:
    void profile(const std::shared_ptr<profile>& p) {
      profile_ = p;
    }
    const flick::profile& profile() const {return *profile_;}
    virtual double absorption_coefficient() const = 0;
    // [per length]

    virtual double scattering_coefficient() const = 0;
    // [per length]

    double scattering_optical_depth(double height, double angle, double distance) {
      return concentration()*scattering_coefficient()*
	profile().column(height,angle,distance);
    }
    double distance(double height, double angle,
		    double scattering_optical_depth) {
      return 0;
    }
  };
  */
  /*
  class composite_material : public material {
    std::vector<std::shared_ptr<fabric>> materials_; 
  public:
    double absorption_coefficient() override {
      double a = 0
      for (size_t i=0; i<materials_.size(); ++i)
	a += materials_[i].absorption_coefficient();
      return a;
    }
    void add(std::shared_ptr<fabric> f) {
      materials_.emplace_back(f);
    }
  };

  {
    composite_material ice;
    auto m = std::make_shared<pure_ice>();
    m->temperature(273);
    m->salinity(0.03);
    ice->add(m)
  }
  class value { // maybe not?
    double reference;
    double value;
  public:
    value(double reference, double v) {
      frequency f( ) NO.
    }
    double dB() {}
  };
  */
}
}

#endif
