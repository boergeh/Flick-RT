#ifndef flick_coating
#define flick_coating

#include "../environment/input_output.hpp"
#include "../numeric/constants.hpp"
#include "../polarization/mueller.hpp"
//#include "../polarization/fresnel.hpp"


namespace flick {
namespace coating {
  
  class angle_generator {
  protected:
    double ui_polar_{0};
    double ui_azimuth_{0};
  };

  class lambert_angle_generator : public angle_generator {
  public:
    double polar_angle() {
      return asin(sqrt(ui_polar_));
    }
    double azimuth_angle() {
      return ui_azimuth_ * 2 * constants::pi;
    }
    double brdf() {
      return 1/(2*constants::pi);
    }
  };
 
  class base {
  public:
    virtual mueller reflection_mueller_matrix() const = 0;
    virtual mueller transmission_mueller_matrix() const = 0;
    virtual double reflectivity() const = 0; 
    virtual double transmissivity() const = 0;
  };

  /*
  class isotropic : public base {
  protected:
    double incidence_angle_{0};
  public:
    void incidence_angle(double ia) {
      incidence_angle_ = ia;
    }
  };
  class basic_lambert : public base {

  };
  */
  
  class grey_lambert : public base {
    double reflectivity_{1};
    double transmissivity_{0};
  public:
    grey_lambert(double reflectivity, double transmissivity)
      : reflectivity_{reflectivity}, transmissivity_{transmissivity}
    {}
    mueller reflection_mueller_matrix() const {
      mueller m;
      m.add(0,0, 1/(2*constants::pi));
      return m;
    }
    mueller transmission_mueller_matrix() const {
      return reflection_mueller_matrix();
    }
    double reflectivity() const {
      return reflectivity_;
    }
    double transmissivity() const {
      return transmissivity_;
    }
  };

  class white_lambert : public grey_lambert {
  public:
    white_lambert() : grey_lambert::grey_lambert(1,0) {}
  };

  class loamy_sand : base {
    double wavelength_{500e-9};
  public:
    void wavelength(double wl) {
      wavelength_ = wl;
    }
    double reflectivity() const {
      // wavelength dependent tbi
      return 0;
    }
    double transmissivity() const {
      // wavelength dependent tbi
      return 0;
    }
  };

  
 
  /*
  // move to transport?
  class fresnel : public isotropic {
    std::complex<double> m_{1,0};
    flick::fresnel f_;
    //double surface_slope_variance_{0}; 
  public:
    fresnel(const std::complex<double>& relative_refractive_index)
      : m_{relative_refractive_index}, f_{m_, 0} {}
    void incidence_angle(double ia) {
      incidence_angle_ = ia;
      f_ = flick::fresnel(m_, ia);
    }
    //void surface_slope_variance(double ssv) {
    //  surface_slope_variance_ = ssv;
    //}
    mueller reflection_mueller_matrix() {
      return reflection_mueller(f_);// *
      //angular_distribution(f_.reflection_angle());
    }
    mueller transmission_mueller_matrix() {
      return transmission_mueller(f_);// *
      //angular_distribution(std::real(f_.transmission_angle()));
    }
    double reflectivity() {
      return f_.R(); 
    }
    double transmissivity() {
      return f_.T();
    }
  };
  */

  /*
  class isotropic_composite : public isotropic {
    std::vector<std::shared_ptr<isotropic>> isotropics_;
    std::vector<double> fractions_;
  public:
  };
  */

}
}

#endif
