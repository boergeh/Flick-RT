#ifndef flick_stokes
#define flick_stokes

#include "../numeric/constants.hpp"

namespace flick {
  class stokes
  // see Wikipedia Stokes parameters
  {
    std::vector<double> s_{1,0,0,0};
  public:
    stokes() = default;
    stokes(double s0, double s1, double s2, double s3)
      : s_{s0,s1,s2,s3} {};
    /*
    stokes(double intensity, double rotation_angle,
	   double eccentricity_angle, double degree_of_polarization) {
      double psi = rotation_angle;
      double xsi = eccentricity_angle;
      double p = degree_of_polarization;
      s_[0] = intensity;
      s_[1] = intensity*p*cos(2*psi)*cos(2*xsi);
      s_[2] = intensity*sin(2*psi)*cos(2*xsi);
      s_[3] = intensity*p*sin(2*xsi);
    }
    */
    stokes(const std::vector<double>& s) : s_{s} {
      ensure(s.size()==4);
    }
    static stokes polarization_ellipse(double intensity, double rotation_angle,
	   double eccentricity_angle, double degree_of_polarization) {
      double psi = rotation_angle;
      double xsi = eccentricity_angle;
      double p = degree_of_polarization;
      double s0 = intensity;
      double s1 = intensity*p*cos(2*psi)*cos(2*xsi);
      double s2 = intensity*sin(2*psi)*cos(2*xsi);
      double s3 = intensity*p*sin(2*xsi);
      return stokes{s0,s1,s2,s3};
    }
    static stokes unpolarized() {
      return stokes{1,0,0,0};
    }
    static stokes s_polarized() {
      return stokes{1,-1,0,0};
    }
    static stokes p_polarized() {
      return stokes{1,1,0,0};
    }
    static stokes lhc_polarized() {
      return stokes{1,0,0,-1};
    }
    static stokes rhc_polarized() {
      return stokes{{1,0,0,1}};
    }
    double operator()(size_t element_number) {
      return s_.at(element_number);
    }
    double I() const {return s_[0];}
    double Q() const {return s_[1];}
    double U() const {return s_[2];}
    double V() const {return s_[3];}
    
    stokes& rotate(double delta_psi) {
      double u = 2 * delta_psi;
      double s1r = cos(u)*s_[1] + sin(u)*s_[2];
      double s2r = -sin(u)*s_[1] + cos(u)*s_[2];
      s_[1] = s1r;
      s_[2] = s2r;
      return *this;
    }
    stokes& scale(double factor) {
      for (size_t i=0; i<s_.size(); ++i)
	s_[i] *= factor;
      return *this;
    }
    double rotation_angle() const {
      double psi = constants::pi/2;
      if (s_[1] != 0)
	psi = 0.5*atan(s_[2]/s_[1]);
      if ((s_[1] < 0 && s_[2] >= 0) || (s_[2] >= 0 && s_[1] < 0))
	return psi + constants::pi/2;
      return psi;
    }
    double eccentricity_angle() const {
      double denominator = sqrt(pow(s_[1],2) + pow(s_[2],2));
      return 0.5*atan(s_[3]/denominator);
    }
    double degree_of_polarization() const {
      double polarized_part = sqrt(pow(s_[1],2) + pow(s_[2],2) + pow(s_[3],2));
      return polarized_part/s_[0];
    }
    friend std::ostream& operator<<(std::ostream &os, const stokes& s) {
      os << "[I Psi Xsi p] " << s.I() << " " << s.rotation_angle()
	 << " " << s.eccentricity_angle()
	 << " " << s. degree_of_polarization() << ", " <<
	"[I Q U V] "<< s.I() << " " << s.Q() << " " << s.U()
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
