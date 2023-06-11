#include "matrix.hpp"

namespace flick {
  begin_test_case(matrix_test) {
    using namespace linalg;
    matrix m1 = {{1,2},{3,4}};
    check_close(det(m1),-2,1e-12,"m1");

    matrix m2 = {{4,3},{1,5}};
    m2 = cofactor(m2);
    check_close(m2[0][0],5,1e-12);
    check_close(m2[0][1],-1,1e-12);
    check_close(m2[1][0],-3,1e-12);
    check_close(m2[1][1],4,1e-12);

    matrix m3 = T(m2);
    check_close(m3[0][0],5,1e-12);
    check_close(m3[1][0],-1,1e-12);

    
    matrix m4 = T(m3);
    check_close(m4[1][0],-3,1e-12);

    matrix m5 = {{4,3},{1,5},{1,3}};
    T(m5);

    matrix m6 = {{4,7},{2,6}};
    m6 = inv(m6);
    check_close(m6[0][0],0.6,1e-12);
    check_close(m6[0][1],-0.7,1e-12);
    check_close(m6[1][0],-0.2,1e-12);
    check_close(m6[1][1],0.4,1e-12);
    
    matrix m7_a = {{1,2,3},{4,5,6}};
    matrix m7_b = {{7,8},{9,10},{11,12}};
    matrix m7 = m7_a * m7_b;
    check_close(m7[0][0],58,1e-12);
    check_close(m7[0][1],64,1e-12);
    check_close(m7[1][0],139,1e-12);
    check_close(m7[1][1],154,1e-12);

    matrix m8 = {{0,1},{1,1},{2,1}};
    std::vector<double> v = {1,2,3};
    std::vector<double> c = solve(m8,v);
    std::cout << c.size() << std::endl;
    check_close(c.at(0),1,1e-12);
    check_close(c.at(1),1,1e-12);
    
      
    
  } end_test_case()  
}
