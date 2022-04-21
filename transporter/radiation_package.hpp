#ifndef flick_radiation_package
#define flick_radiation_package

#include "../numeric/direction_generator.hpp"

namespace flick {
  using namespace constants;
  class stokes
  // see Wikipedia Stokes parameters
  {
    double I_{1}; 
    double psi_{0};
    double xsi_{pi/2};
    double p_{0};
  public:
    stokes() = default;
    stokes(double intensity, double psi_angle,
	   double xsi_angle, double degree_of_polarization) 
      : I_{intensity}, psi_{psi_angle}, xsi_{xsi_angle},
	p_{degree_of_polarization} {}
    stokes(const std::vector<double>& s) {
      double s0 = s.at(0);
      double s1 = s.at(1);
      double s2 = s.at(2);
      double s3 = s.at(3);
      I_ = s0;
      if (s0 > 0)
	p_ = sqrt(s1*s1 + s2*s2 + s3*s3)/s0;
      else
	p_ = 0;
      if (s1 > 0)
	psi_ = 0.5*atan(s2/s1);
      else
	psi_ = 0;
      double denominator = sqrt(s1*s1+s2*s2);
      if (denominator > 0)
	xsi_ = 0.5*atan(s3/denominator);
      else
       xsi_ = 0;
    }
    double I() const {return I_;}
    double Q() const {return I_*p_*cos(2*psi_)*cos(2*xsi_);}
    double U() const {return I_*p_*sin(2*psi_)*cos(2*xsi_);}
    double V() const {return I_*p_*sin(2*xsi_);}
    std::vector<double> vector() const {return std::vector{I(),Q(),U(),V()};}
    stokes& rotate_ellipse(double delta_psi) {
      psi_ += delta_psi;
      return *this;
    }
    stokes& I(double i) {I_ = i; return *this;}
    stokes& psi_angle(double pa) {
      psi_ = pa;
      return *this;
    }
    double psi_angle() const {return psi_;}
    double xsi_angle() const {return xsi_;}
    double degree_of_polarization() const {return p_;}
    friend std::ostream& operator<<(std::ostream &os, const stokes& s) {
      os << s.I_ << " " << s.psi_ << " " << s.xsi_ << " " << s.p_;
      return os;
    }
  };

  class radiation_package
  // First stokes parameter is intensity weight. All radiation
  // variables may change during transport. All stokes
  // parameters have units of per solid angle.
  {
    pose pose_;
    double wavelength_{500e-9};
    stokes stokes_;
    double traveling_length_{0};
    size_t scattering_events_{0};
    //bool do_not_scatter_?
  public:
    radiation_package() = default;
    radiation_package(double wavelength, const stokes& s)
      : wavelength_{wavelength}, stokes_{s} {}
    void scale_intensity(double factor) {
      stokes_.I(stokes_.I()*factor);
    }
    void move_to(const vector& position) {
      pose_.move_to(position);
    }
    void rotate_to(const unit_vector& direction) {//?? need more directional information
      pose_.rotate_to(direction);
    }
    auto pose(){
      return pose_;
    }
    auto wavelength() {
      return wavelength_;
    }
    friend std::ostream& operator<<(std::ostream &os, const radiation_package& rp) {
      os << rp.pose_ << " " << rp.wavelength_ << " " << rp.stokes_ << " "
	 << rp.traveling_length_ << " " << rp.scattering_events_;
      return os;
    }
  };
}
#endif
