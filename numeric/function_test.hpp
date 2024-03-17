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
    check_throw(fk.value(0));
    pe_function fl(2);
    check_close(fl.value(),2);
    pp_function fm{{1, 2, 3, 4},{1,4,9,16}};
    fm.add_extrapolation_points(1);
    check_close(fm.value(0.01),1,1e-9_pct);    
    fm.scale_x(2);
    check_close(fm.value(2),1);
    fm.scale_y(2);
    check_close(fm.value(2),2);
    pp_function fmb = significant_digits(fm, 10, 1);
    check_close(fmb.value(8),30);
    check_throw(fmb.add_extrapolation_points(0));
    
    function<piecewise_linear> fn{{-2,-1,5},{-2,-1,5}};
    check_close(fn.value(-3),-3);
    check_close(fn.derivative(0),1);
    check_close(fn.integral(-4,6),10);
    check_close(*fn.integral_limit_b(-4,10),6);
    check_close(*fn.integral_limit_b(6,-10),-4);
    fn.add_extrapolation_points(0);
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
    check_close(fa.integral(),1);
    pp_function fb{{0.1,1,2,2.5,3,4},{0,0,1,1,1,0}};
    check_small(fb.value(1));   
    check_close(fb.value(3),1);   
    check_close(fb.value(2),1);   
    check_small(fb.value(3.5));
    check_close(fb.integral(),1);
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
}
