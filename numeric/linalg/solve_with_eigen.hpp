#ifndef flick_linalg_solve_with_eigen
#define flick_linalg_solve_with_eigen

#include "matrix.hpp"
#include "Eigen/Dense"

namespace flick {
  namespace linalg {  
    std::vector<double> solve(const matrix& m, const std::vector<double>& v)
    // Solves set of linear equations using the Eigen library
    {
      int rows = m.size();
      int cols = m.at(0).size();
      Eigen::MatrixXd me(rows,cols);
      for (size_t i=0; i<rows; i++)
	for (size_t j=0; j<cols; j++)
	  me(i,j) = m[i][j];
      const Eigen::VectorXd b = Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(v.data(),v.size());
      Eigen::VectorXd a = me.colPivHouseholderQr().solve(b);
      return std::vector<double>(a.data(),a.data()+a.size());
    }
  }
}

#endif
