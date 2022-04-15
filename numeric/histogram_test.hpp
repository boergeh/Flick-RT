#include "histogram.hpp"
namespace flick {
  using namespace constants;
  begin_test_case(histogram_test) {  
    equal_bins x_bins{-1, 1, 2};
    equal_bins y_bins{-1, 1, 3};
    histogram h(x_bins,y_bins);
    histogram h1d(x_bins);
    h.add(-0.01,0,pi);
    h.add(h);
    check_small(h.bin_value(0,1)-2*pi,1e-15);
    write(h,"histogram_view.txt");
  } end_test_case()
}
