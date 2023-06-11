#ifndef flick_linalg_matrix
#define flick_linalg_matrix

#include <iostream>
#include <vector>

namespace flick {
  namespace linalg
  // Most of these methods are made with help from ChatGPT-2023
  {
    using matrix = std::vector<std::vector<double>>;

    std::ostream& operator<<(std::ostream& os, const matrix& m) {
      os << "\n";
      for (size_t i=0; i < m.size(); i++) {
	for (size_t j=0; j < m.at(0).size(); j++) {
	  os << m[i][j] << " ";
	}
	os << "\n";
      }
      return os;
    }

    matrix operator*(double k, matrix m) {
      for (size_t i=0; i < m.size(); i++) {
	for (size_t j=0; j < m[0].size(); j++) {
	  m[i][j] = k * m[i][j];
	}
      }
      return m;
    }
    
    matrix operator*(const matrix& m, double k) {
      return k * m;
    }

    matrix operator*(const matrix& m1, const matrix& m2)
    // Matrix multiplication
    {
      size_t rows1 = m1.size();
      size_t cols1 = m1[0].size();
      size_t rows2 = m2.size();
      size_t cols2 = m2[0].size();
      if (cols1 != rows2) {
	std::cerr <<
	  "Error: Invalid matrix dimensions for multiplication!"
		  << std::endl;
	return matrix();
      }
      matrix result(rows1, std::vector<double>(cols2, 0.0));
      for (size_t i = 0; i < rows1; i++) {
	for (size_t j = 0; j < cols2; j++) {
	  for (size_t k = 0; k < cols1; k++) {
	    result[i][j] += m1[i][k] * m2[k][j];
	  }
	}
      }
      return result;
    }

    double det(const matrix& m)
    // Determinant
    {
      size_t n = m.size();
      if (n == 1) {
	return m[0][0];
      }
      double d = 0.0;
      int sign = 1;
      for (size_t col = 0; col < n; col++) {
	matrix submatrix(n - 1, std::vector<double>(n - 1));
	size_t sub_i = 0;
	for (size_t row = 1; row < n; row++) {
	  size_t sub_j = 0;
	  for (size_t subcol = 0; subcol < n; subcol++) {
	    if (subcol == col) continue;
	    submatrix[sub_i][sub_j] = m[row][subcol];
	    sub_j++;
	  }
	  sub_i++;
	}
	double subDet = det(submatrix);
	d += sign * m[0][col] * subDet;
	sign *= -1;
      }
      return d;
    }

    matrix cofactor(const matrix& m) {
      size_t n = m.size();
      matrix cofactorMatrix(n, std::vector<double>(n));

      for (size_t i = 0; i < n; i++) {
	for (size_t j = 0; j < n; j++) {
	  matrix submatrix(n - 1, std::vector<double>(n - 1));
	  size_t sub_i = 0;
	  for (size_t row = 0; row < n; row++) {
	    if (row == i) continue;
	    size_t sub_j = 0;
	    for (size_t col = 0; col < n; col++) {
	      if (col == j) continue;
	      submatrix[sub_i][sub_j] = m[row][col];
	      sub_j++;
	    }
	    sub_i++;
	  }
	  double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
	  double cofactor = sign * det(submatrix);
	  cofactorMatrix[i][j] = cofactor;
	}
      } 
      return cofactorMatrix;
    }

    matrix T(const matrix& m)
    // Transpose
    {
      size_t rows = m.size();
      size_t cols = m[0].size();
      matrix t(cols, std::vector<double>(rows));
      for (size_t i = 0; i < rows; i++) {
	for (size_t j = 0; j < cols; j++) {
	  t[j][i] = m[i][j];
	}
      }
      return t;
    }

    matrix inv(const matrix& m)
    // Matrix inversion
    {
      return 1/det(m)*T(cofactor(m));
    }
  }
}

#endif
