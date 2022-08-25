#ifndef flick_stokes
#define flick_stokes

#include "../numeric/constants.hpp"

namespace flick {
  class stokes
  // see Wikipedia Stokes parameters
  {
    double I_{1}; //should store log(I)?
    double psi_{0};
    double xsi_{constants::pi/2};
    double p_{0};
  public:
    stokes() = default;
    stokes(double intensity, double rotation_angle,
	   double eccentricity_angle, double degree_of_polarization) 
      : I_{intensity}, psi_{rotation_angle}, xsi_{eccentricity_angle},
	p_{degree_of_polarization} {
    }
    stokes(const std::vector<double>& s) : p_{0}, psi_{0}, xsi_{0} {
      ensure(s.size()==4);
      double polarized_part = sqrt(pow(s[1],2) + pow(s[2],2) + pow(s[3],2));
      I_ = s[0];
      if (s[0] > 0)
	p_ = polarized_part/s[0];
      if (s[1] > 0)
	psi_ = 0.5*atan(s[2]/s[1]);
      double denominator = sqrt(pow(s[1],2) + pow(s[2],2));
      if (denominator > 0)
	xsi_ = 0.5*atan(s[3]/denominator);
    }
    double I() const {return I_;}
    double Q() const {return I_*p_*cos(2*psi_)*cos(2*xsi_);}
    double U() const {return I_*p_*sin(2*psi_)*cos(2*xsi_);}
    double V() const {return I_*p_*sin(2*xsi_);}
    stokes& rotate(double delta_psi) {
      psi_ += delta_psi;
      return *this;
    }
    stokes& scale(double factor) {
      I_ *= factor;
      return *this;
    }
    double rotation_angle() const {return psi_;}
    double eccentricity_angle() const {return xsi_;}
    double degree_of_polarization() const {return p_;}
    friend std::ostream& operator<<(std::ostream &os, const stokes& s) {
      os << "[I Psi Xsi p]: " << s.I_ << " " << s.psi_ << " " << s.xsi_
	 << " " << s.p_ << ", " <<
	"[I Q U V]: "<< s.I() << " " << s.Q() << " " << s.U()
	 << " " << s.V() <<'\n';
      return os;
    }
  private:
    void ensure(bool b) {
      if (!b)
	throw std::invalid_argument("stokes");
    }
  };
}

#endif
