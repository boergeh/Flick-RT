#ifndef flick_mueller
#define flick_mueller

namespace flick {
  class mueller
  // Holds non-zero Mueller matrix element values
  {
    struct element {
      size_t row;
      size_t col;
      double value;
    };
    std::vector<element> elements_;
  public:
     mueller& add(size_t row, size_t col, double value) {
      elements_.push_back(element{row,col,value});
      return *this;
    }
    size_t size() const {
      return elements_.size();
    }
    const element& operator()(size_t n) const {
      return elements_.at(n);
    }
  };
  
  std::ostream& operator<<(std::ostream &os, const mueller& m) {
    size_t n = 0;
    os << '\n';
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
	if (n < m.size() && m(n).row == i && m(n).col == j) {
	  os << m(n).value << " ";
	  n++;
	} else {
	  os << 0 << " ";
	}
      }
      os << '\n';
    }
    return os;
  }
}

#endif
