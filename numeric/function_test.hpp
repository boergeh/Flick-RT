#include "function.hpp"

namespace flick {
  begin_test_case(function_test_A) {
    function<piecewise_exponential> fa{{-1,1},{1,1}};
    check_small(fa.derivative(-3));
    check_close(fa.value(-3), 1);
    check_close(fa.integral(-3,-2), 1);
    check_close(fa.integral(-3,3), 6);
    double integral_value = fa.integral(-3,-2);
    std::optional<double> limit = fa.integral_limit_b(-3,integral_value);
    check_close(*limit, -2);
    integral_value = fa.integral(-2, 2);   
    limit = fa.integral_limit_b(-2,integral_value);
    check_close(*limit, 2);

    function<piecewise_exponential> fb{{-1,0,1,1.1},
				       {exp(-2*1),exp(0),exp(2*1),exp(2*1.1)}};
    check_close(fb.value(-1.01), exp(-1.01*2));
    check_close(fb.derivative(-1.01), exp(-1.01*2)*2);

    function<piecewise_exponential> fc{{-1,0,1,1.1},
				       {exp(-1),exp(0),exp(1),exp(1.1)}};
    check_close(fc.value(-3), exp(-3));
    check_close(fc.value(3), exp(3));
    check_close(fc.derivative(3), exp(3));
    double area_c4 = exp(3)-exp(-3);
    check_close(fc.integral(-3,3), area_c4);
    check_close(fc.integral(3,-3), -area_c4);
    check_close(fc.integral(-1,1), exp(1)-exp(-1));
    check_close(fc.integral(-1,0.5), exp(0.5)-exp(-1));    
    integral_value = fc.integral(-3,-2);
    limit = fc.integral_limit_b(-3,integral_value);
    check_close(*limit, -2);
    integral_value = fc.integral(-1,1);
    limit = fc.integral_limit_b(-1,integral_value);
    check_close(*limit, 1);
    integral_value = fc.integral(-2,2);
    limit = fc.integral_limit_b(-2,integral_value);
    check_close(*limit, 2);
    
    limit = fc.integral_limit_b(2,-integral_value);
    check_close(*limit, -2); 
    integral_value = fc.integral(2,3);
    limit = fc.integral_limit_b(2,integral_value);
    check_close(*limit, 3);
    limit = fc.integral_limit_b(3,-integral_value);
    check_close(*limit, 2);
    
    function<piecewise_exponential> fd{{-1,0,1,2},
				       {exp(1),exp(0),exp(-1),exp(-2)}};
    check_close(fc.integral(-1,1.9), exp(1.9)-exp(-1));
    integral_value = 10;
    limit = fd.integral_limit_b(-2,integral_value);
    check(!limit.has_value());
    
    function<piecewise_power> fe{{1,2},{1,1}};
    check_small(fe.derivative(3));
    check_close(fe.value(3), 1);
    check_close(fe.integral(3,7), 4);
    check_close(fe.integral(0.5,3), 2.5);

    function<piecewise_power> ff{{exp(-1),1},{1, exp(-1)}};
    double area = exp(-1)*log(3/0.5);
    check_close(ff.integral(0.5,3), area);
    limit = ff.integral_limit_b(2,1000);
    check(limit.has_value());
    integral_value = ff.integral(0.1,10);
    limit = ff.integral_limit_b(0.1,integral_value);
    check_close(*limit, 10);

    function<piecewise_power> fg{{0.1,1,7},{pow(0.1,2),pow(1,2),pow(7,2)}};
    check_small(fg.integral(1,1));
    check_close(fg.value(1.1), pow(1.1,2));
    check_close(fg.derivative(7), 14);
    check_close(fg.integral(1,2), 1./3*(pow(2,3)-1));
    integral_value = fg.integral(0.1,10);
    limit = fg.integral_limit_b(0.1,integral_value);
    check_close(*limit, 10);
    limit = fg.integral_limit_b(10,-integral_value);
    check_close(*limit, 0.1, 1e-8_pct);

    function<piecewise_power> fh{{0.1,0.5,1,7,8,12},{1,0.1,2,3,1,6}};
    integral_value = fh.integral(0.01,10);
    limit = fh.integral_limit_b(0.01,integral_value);
    check_close(*limit, 10);
    limit = fh.integral_limit_b(10,-integral_value);
    check_close(*limit, 0.01, 1e-9_pct);
    function<piecewise_linear> fh_l{{0.1,0.5,1,7,8,12},{1,0.1,2,3,1,6}};
    check_small(*fh_l.integral_limit_b(10,-integral_value), 0.2);
    
    pp_function fi{{1,2},{1,1}};
    pe_function fj{{1,2},{1,1}};
    pe_function fk;
    check_throw(pp_function({-1,2},{1,1}));
    pe_function fl(2);
    check_close(fl.value(),2);
    
    pp_function fm{{1, 2, 3, 4},{1,4,9,16}};
    fm = fm.constant_extrapolation();
    check_close(fm.value(0.01),1,1e-9_pct);    
    fm.scale_x(2);
    check_close(fm.value(2),1);
    fm.scale_y(2);
    check_close(fm.value(2),2);
    pp_function fmb = significant_digits(pp_function{{1, 2, 3, 4},{1,4,9,16}}, 10, 1);
    check_close(fmb.value(4),20);
 
    function<piecewise_linear> fn{{-2,-1,5},{-2,-1,5}};
    check_close(fn.value(-3),-3);
    check_close(fn.derivative(0),1);
    check_close(fn.integral(-4,6),10);
    check_close(*fn.integral_limit_b(-4,10),6);
    check_close(*fn.integral_limit_b(6,-10),-4);
    fn = fn.zero_extrapolation();
    //fn.add_extrapolation_points(0);
    check_small(fn.value(-5));
       
    function<piecewise_linear> fo{{-2,-1,5},{0,0,0}};
    check_small(fo.value(-3));
    check_small(fo.derivative(0));
    check_small(fo.integral(-4,6));
    check(fo.integral_limit_b(-4,10).has_value()==false);

    function<piecewise_linear> fp{{-1,0,1},{0,1,2}};
    auto acc = fp.accumulation();
    check_close(acc.back(),2);
    auto icd = inverted_cumulative_distribution(fp);
    check_close(icd.value(1),1);
    auto is = importance_sampled(fp,10);
    check_close(is.integral(),fp.integral());

    pp_function ppf{{400e-9,400.1e-9},{1,1e6}};
    check(!std::isnan(ppf.value(400.05e-9)));
    check(!std::isnan(ppf.integral(399e-9,401e-9)));
    check(!std::isnan(ppf.derivative(399e-9)));
    check(ppf.integral_limit_b(399e-9,1).has_value());

    pe_function pef{{400e-9,400.1e-9},{1,1e6}};
    check(!std::isnan(pef.value(400.05e-9)));
    check(!std::isnan(pef.integral(399e-9,401e-9)));
    check(!std::isnan(pef.derivative(399e-9)));
    check(pef.integral_limit_b(399e-9,1).has_value());

    pe_function peg{{1},{1}};
    check_small(peg.integral());
    check_close(peg.integral(-1,1),2);
  } end_test_case()
   
  begin_test_case(function_test_B) {
    pl_function fa{{0,1,2},{0,1,0}};
    check_close(fa.integral(),1);
    pl_function fb{{-2,-1,0},{0,-1,0}};
    check_close(fb.integral(),-1);   
  } end_test_case()

  begin_test_case(function_test_C) {
    pl_function fa{{-1,0,2},{-1,1,2}};
    pl_function fb{{0,2},{1,2}};
    std::vector<double> x = {-1.1, 0, 0.8};
    pl_function fc = integral_conservative_add(fa,fb,x);
    check_close(fa.integral()+fb.integral(),fc.integral());   
  } end_test_case()
  
  begin_test_case(function_test_D) {
    pe_function fa{{-2,-1,0,1,2},{0,0,1,1,0}};
    check_small(fa.value(-3));   
    check_small(fa.value(-2.5));   
    check_small(fa.value(-0.9));   
    check_close(fa.value(1),1);   
    check_close(fa.value(1),1);   
    check_small(fa.value(1.1));
    check_close(fa.integral(),1,0.5_pct);
    pp_function fb{{0.1,1,2,2.5,3,4},{0,0,1,1,1,0}};
    check_small(fb.value(1));   
    check_close(fb.value(3),1);   
    check_close(fb.value(2),1);   
    check_small(fb.value(3.5));
    check_close(fb.integral(),1,0.5_pct);
  } end_test_case()
   
  begin_test_case(function_test_E) {
    std::istringstream s("/* test stream */ 1 1 2 2 3 3");
    pl_function f;
    s >> f;
    check_close(f.integral(0,10),50);
  } end_test_case()
  
  begin_test_case(function_test_F) {
    // Integration in backward direction
    std::vector<double> x = {-3, -2, -1, 0};
    std::vector<double> y = {0.1, 0.1, 1, 1};
    pl_function f(x,y);
    check_close(f.integral(0,-1),-1);
    pe_function g(x,y);
    check_close(g.integral(0,-1),-1);
  } end_test_case()
  
  begin_test_case(function_test_G) {
    // Curvature
    std::vector<double> x = {-3, -2, -1, 0, 1};
    std::vector<double> y = {exp(-3), exp(-2), exp(-1), exp(0), exp(1)};
    pe_function f{x,y};
    pl_function c = derivative(to_exponential(derivative(f))); 
    check_close(c.value(-3),exp(-3));
    pl_function pdf = absolute(c).normalize();
    check_close(pdf.integral(),1);
    pl_function cdf = accumulate(pdf);
    check_small(cdf.y().front());
    check_close(cdf.y().back(),1);
    cdf = remove_non_increasing_values(cdf);
    pl_function quantile = invert(cdf);
    check_close(quantile.value(0),-3);     
  } end_test_case()
  
  begin_test_case(function_test_H) {
    // Extrapolation
    std::vector<double> x = {1, 1e9};
    std::vector<double> y = {1, 1};
    pp_function f{x,y};
    pp_function g = f.zero_extrapolation();
    check_close(g.integral(1e-19,1e19),1e9-1);
    pe_function f2{x,y};
    pe_function g2 = f2.zero_extrapolation();
    check_close(g2.integral(1e-19,1e19),1e9-1);
    pe_function f3{{-1,1},{0,0}};
    pe_function g3 = f3.zero_extrapolation();
    check_small(g3.integral(1e-19,1e19));
    pe_function f4{{-1,1},{1,1}};
    pe_function g4 = f4.constant_extrapolation();
    check_close(g4.integral(-1e19,1e19),2*1e19);
  } end_test_case()
  
  begin_test_case(function_test_I) {
    // Small steps
    double epsilon = 1e-6;
    std::vector<double> x = {-100, -epsilon};
    std::vector<double> y = {1e30, 1e30};
    pl_function f{x,y};
    f = f.zero_extrapolation();
    check_close(f.integral(-101,1),f.integral(-100,-epsilon),1e-10);
    std::vector<double> x2 = {-100, 0, 1e-3*epsilon, epsilon};
    std::vector<double> y2 = {1e30, 1e30, 1e30, 1e30};
    pe_function f2{x2,y2};
    f2 = f2.zero_extrapolation();
    check_close(f2.integral(-101,1),f2.integral(-100,epsilon));
    std::vector<double> x3 = {epsilon, 2*epsilon, 0.1};
    std::vector<double> y3 = {1e30, 1e30, 1e30};
    pp_function f3{x3,y3};
    f3 = f3.zero_extrapolation();
    double x_low = 0.1*epsilon;
    check_close(f3.integral(x_low,1),f3.integral(x3.front(),x3.back()),1e-10);
    std::vector<double> x4 = {1, 2};
    std::vector<double> y4 = {1e308, 1e308};
    pl_function f4{x4,y4};
    f4 = f4.zero_extrapolation();
    check_close(f4.value(-1),f4.value(-1e9));
    double dx = 9.99999999251599547e-07;
    double y5 = 0.02;
    pl_function f5{{0, dx},{y5,y5}};
    f5 = f5.zero_extrapolation();
    check_close(f5.integral(0,dx),y5*dx);   
  } end_test_case()
  
  begin_test_case(function_test_J) {
    double y1 = 0.0208156402400279592;
    double y2 = 0.0208156402400238166;
    double x1 = -200;
    double x2 = 0;
    double a = -100;
    double b = -99.9;
    pe_function f{{x1, x2},{y1,y2}};
    check_close(f.integral(a,b),(f.value(a)+f.value(b))/2*(b-a));   
  } end_test_case()
}
