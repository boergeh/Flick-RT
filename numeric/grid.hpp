#ifndef flick_grid
#define flick_grid

#include <vector>

namespace flick {
  class grid_4d {
    size_t nd_ = 4;
  public:
    std::vector<std::vector<double>> x;
    std::vector<std::vector<std::vector<std::vector<double>>>> f;
    grid_4d() = default;
    grid_4d(int ni, int nj, int nk, int nl) {
      x.resize(nd_);
      x[0].resize(ni);
      x[1].resize(nj);
      x[2].resize(nk);
      x[3].resize(nl);
      f.resize(ni);
      for (int i = 0; i < ni; i++) {
	f[i].resize(nj);
	for (int j = 0; j < nj; j++) {
	  f[i][j].resize(nk);
	  for (int k = 0; k < nk; k++) {
	    f[i][j][k].resize(nl);
	  }
	}
      }
    }
  private:
    friend std::ostream& operator<<(std::ostream &os, const grid_4d& f) {
      for (size_t i=0; i<f.x.size(); i++) {
	os << f.x[i].size() << " ";
      }
      os << "\n";
      for (size_t i=0; i<f.x.size(); i++) {
	os << f.x[i] << "\n";
      }
      for (int i = 0; i < f.x[0].size(); i++) {
	for (int j = 0; j < f.x[1].size(); j++) {
	  for (int k = 0; k < f.x[2].size(); k++) {
	    for (int l = 0; l < f.x[3].size(); l++) {
	      os << f.f[i][j][k][l] << " ";
	    }
	  }
	}
      }
      os << "\n";
      return os;
    }
    friend std::istream& operator>>(std::istream &is, grid_4d& f) {
      std::vector<int> n(f.nd_);
      for (size_t i=0; i<f.nd_; i++) {
	is >> n[i];
      }
      f = grid_4d(n[0],n[1],n[2],n[3]);
      for (size_t i=0; i<f.nd_; i++) {
	for (size_t j=0; j<n[i]; j++) {
	  is >> f.x[i][j];
	}
      }
      for (int i = 0; i < n[0]; i++) {
	for (int j = 0; j < n[1]; j++) {
	  for (int k = 0; k < n[2]; k++) {
	    for (int l = 0; l < n[3]; l++) {
	      is >> f.f[i][j][k][l];
	    }
	  }
	}
      }
      return is;
    }
  };
}

#endif
