#include "spherical_harmonics.hpp"

namespace flick {
  begin_test_case(spherical_harmonics_test) {
    using Y = spherical_harmonics;
    using complex = std::complex<double>;
    double pi = std::numbers::pi;
    complex i(0,1);
    direction d{0.1,0.1};
    size_t l_max = 1;  
    complex Y_00(1./2*sqrt(1/pi),0);
    complex Y_1minus1 = 1./2*sqrt(3/(2*pi))*sin(d.theta)*exp(-d.phi*i);
    check_close(std::real(Y(d,l_max).lm(0,0)),std::real(Y_00));
    check_close(std::real(Y(d,l_max+1).lm(0,0)),std::real(Y_00));
    check_close(std::real(Y(d,l_max).lm(1,-1)),std::real(Y_1minus1));
    check_close(std::imag(Y(d,l_max+1000).lm(1,-1)),std::imag(Y_1minus1));

  } end_test_case() 
}
