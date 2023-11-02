#ifndef flick_linalg_matrix
#define flick_linalg_matrix

#include <iostream>
#include <vector>

namespace flick {
  namespace linalg
  {
    using matrix = std::vector<std::vector<double>>;
    using vector = std::vector<double>;
    
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

    std::istream& operator>>(std::istream &is, matrix& m) {
      m.clear();
      double x;
      std::vector<double> row;
      std::string line;
      while (std::getline(is, line)) {	  
	std::stringstream ss(line);
	while (ss >> x) {
	  row.push_back(x);
	}
	if (row.size() > 0)
	  m.push_back(row);
	row.clear();
      }
      return is;
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
 
    double det(matrix m)
    // Determinant
    {
      int N = static_cast<int>(m.size());
      double d = 1;
      for (int i = 0; i < N; ++i) {	
        double pivotElement = m[i][i];
        int pivotRow = i;
        for (int row = i + 1; row < N; ++row) {
	  if (std::abs(m[row][i]) > std::abs(pivotElement)) {
	    pivotElement = m[row][i];
	    pivotRow = row;
	  }
        }
        if (pivotRow != i) {
	  m[i].swap(m[pivotRow]);
	  d *= -1.0;
        }
        d *= pivotElement;
        for (int row = i + 1; row < N; ++row) {
	  for (int col = i + 1; col < N; ++col) {
	    m[row][col] -= m[row][i] * m[i][col] / pivotElement;
	  }
        }
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

    matrix t(const matrix& m)
    // Transpose
    {
      size_t rows = m.size();
      size_t cols = m[0].size();
      matrix tm(cols, std::vector<double>(rows));
      for (size_t i = 0; i < rows; i++) {
	for (size_t j = 0; j < cols; j++) {
	  tm[j][i] = m[i][j];
	}
      }
      return tm;
    }

    matrix inv(const matrix& m)
    // Matrix inversion
    {
      return 1/det(m)*t(cofactor(m));
    }

    std::vector<double> column(size_t n, const matrix& m) {
      std::vector<double> c(m.size());
      for (size_t i=0; i<m.size(); i++)
	c[i] = m[i][n];
      return c;
    }

    matrix remove_column(size_t n, matrix m) {
      for (size_t i=0; i<m.size(); i++) {
	m[i].erase(m[i].begin() + n);
      }
      return m;
    }
    
    std::vector<double> solve_small(const matrix& m, const std::vector<double>& v)
    // Solves set of linear equations for small matrixes
    {
      auto col = t(matrix{v});
      return column(0, inv(t(m)*m) * (t(m)*col));
    }
  }
}

#endif
