#ifndef flick_coating
#define flick_coating

#include "../environment/input_output.hpp"
#include "../numeric/constants.hpp"
#include "../polarization/fresnel.hpp"
//#include "../polarization/isotropic_mueller.hpp"

namespace flick {
namespace coating {
  class base
  // Note that methods returning reflection and transmission angles
  // are not necessarily giving correct energy distributions. They
  // should rather be used as proposed angles and intensities should
  // be scaled according to transmissivity and reflectivity.
  {
  protected:
    double wavelength_{500e-9};
    double rand_{1};
  public:
    void wavelength(double wl) {
      wavelength_ = wl;
    }
    void random_number(double rn) {
      rand_ = rn;
    }
    virtual mueller reflection_mueller_matrix()=0;
    virtual mueller transmission_mueller_matrix()=0;
    virtual double reflectivity()=0; 
    virtual double transmissivity()=0;
    virtual double reflection_polar_angle()=0;
    virtual double reflection_azimuth_angle()=0;
    virtual double transmission_polar_angle()=0;
    virtual double transmission_azimuth_angle()=0;
  };
  
  class isotropic : public base {
  protected:
    double incidence_angle_{0};
  public:
    void incidence_angle(double ia) {
      incidence_angle_ = ia;
    }
  };

  class basic_lambert : public isotropic {
  public:
    mueller reflection_mueller_matrix() {
      mueller m;
      m.add(0,0, 1/(2*constants::pi));
      return m;
    }
    mueller transmission_mueller_matrix() {
      return reflection_mueller_matrix();
    }
    double reflection_polar_angle() {
      return asin(sqrt(rand_));
    }
    double reflection_azimuth_angle() {
      return rand_ * 2 * constants::pi;
    }
    double transmission_polar_angle() {
      return reflection_polar_angle();
    }
    double transmission_azimuth_angle() {
      return reflection_azimuth_angle();
    }
  };
  
  class grey_lambert : public basic_lambert {
    double reflectivity_{1};
    double transmissivity_{0};
  public:
    grey_lambert(double reflectivity, double transmissivity)
      : reflectivity_{reflectivity}, transmissivity_{transmissivity}
    {}
    double reflectivity() {
      return reflectivity_; 
    }
    double transmissivity() {
      return transmissivity_;
    }
  };
  
  class white_lambert : public grey_lambert {
  public:
    white_lambert() : grey_lambert::grey_lambert(1,0) {}
  };

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
    double reflection_polar_angle() {
      return f_.reflection_angle();
    }
    double reflection_azimuth_angle() {
      return 0;
    }
    double transmission_polar_angle() {
      return std::real(f_.transmission_angle());
    }
    double transmission_azimuth_angle() {
      return 0;
    }
    //private:
    /*
    double angular_distribution(double angle_0) {
      double width = 1e-9;
      double delta_mu = fabs(cos(angle_0)-cos(incidence_angle_));
      std::cout << delta_mu;
      if (delta_mu > width)
	return 0;
      return 1/(2*constants::pi*width);
    }
    */
  };
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
