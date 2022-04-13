#include "function.hpp"

namespace flick {
  begin_test_case(function_test) {
    double epsilon = 1e-12;
    function<piecewise_exponential> fa{{-1,1},{1,1}};
    check_small(fa.derivative(-3), epsilon, "a1");
    check_close(fa.value(-3), 1, epsilon, "a2");
    check_close(fa.integral(-3,-2), 1, epsilon, "a3");
    check_close(fa.integral(-3,3), 6, epsilon, "a4");
    double integral_value = fa.integral(-3,-2);
    std::optional<double> limit = fa.integral_limit_b(-3,integral_value);
    check_close(*limit, -2, epsilon, "a5");
    integral_value = fa.integral(-2, 2);   
    limit = fa.integral_limit_b(-2,integral_value);
    check_close(*limit, 2, epsilon, "a6");

    function<piecewise_exponential> fb{{-1,0,1,1.1},
				       {exp(-2*1),exp(0),exp(2*1),exp(2*1.1)}};
    check_close(fb.value(-1.01), exp(-1.01*2), epsilon, "b1");
    check_close(fb.derivative(-1.01), exp(-1.01*2)*2, epsilon, "c2");

    function<piecewise_exponential> fc{{-1,0,1,1.1},
				       {exp(-1),exp(0),exp(1),exp(1.1)}};
    check_close(fc.value(-3), exp(-3), epsilon, "c1");
    check_close(fc.value(3), exp(3), epsilon, "c2");
    check_close(fc.derivative(3), exp(3), epsilon, "c3");
    double area_c4 = exp(3)-exp(-3);
    check_close(fc.integral(-3,3), area_c4, epsilon, "c4");
    check_close(fc.integral(3,-3), -area_c4, epsilon, "c4b");
    check_close(fc.integral(-1,1), exp(1)-exp(-1), epsilon, "c5");
    check_close(fc.integral(-1,0.5), exp(0.5)-exp(-1), epsilon, "c6");    
    integral_value = fc.integral(-3,-2);
    limit = fc.integral_limit_b(-3,integral_value);
    check_close(*limit, -2, epsilon, "c7");
    integral_value = fc.integral(-1,1);
    limit = fc.integral_limit_b(-1,integral_value);
    check_close(*limit, 1, epsilon, "c8");
    integral_value = fc.integral(-2,2);
    limit = fc.integral_limit_b(-2,integral_value);
    check_close(*limit, 2, epsilon, "c9");
    
    limit = fc.integral_limit_b(2,-integral_value);
    check_close(*limit, -2, epsilon, "c9b"); 
    integral_value = fc.integral(2,3);
    limit = fc.integral_limit_b(2,integral_value);
    check_close(*limit, 3, epsilon, "c10");
    limit = fc.integral_limit_b(3,-integral_value);
    check_close(*limit, 2, epsilon, "c10b");
    
    function<piecewise_exponential> fd{{-1,0,1,2},
				       {exp(1),exp(0),exp(-1),exp(-2)}};
    check_close(fc.integral(-1,1.9), exp(1.9)-exp(-1), epsilon, "d2");
    integral_value = 10;
    limit = fd.integral_limit_b(-2,integral_value);
    check(!limit.has_value(), "d3");
    
    function<piecewise_power> fe{{1,2},{1,1}};
    check_small(fe.derivative(-3), epsilon, "e1");
    check_close(fe.value(-3), 1, epsilon, "e2");
    check_close(fe.integral(3,7), 4, epsilon, "e3");
    check_close(fe.integral(0.5,3), 2.5, epsilon, "e4");

    function<piecewise_power> ff{{exp(-1),1},{1, exp(-1)}};
    double area = exp(-1)*log(3/0.5);
    check_close(ff.integral(0.5,3), area, epsilon, "f5");
    limit = ff.integral_limit_b(2,1000);
    check(limit.has_value(),"f6");
    integral_value = ff.integral(0.1,10);
    limit = ff.integral_limit_b(0.1,integral_value);
    check_close(*limit, 10, epsilon, "f7");

    function<piecewise_power> fg{{0.1,1,7},{pow(0.1,2),pow(1,2),pow(7,2)}};
    check_small(fg.integral(1,1),epsilon,"g0");
    check_close(fg.value(1.1), pow(1.1,2), epsilon, "g1");
    check_close(fg.derivative(7), 7./2, epsilon, "g2");
    check_close(fg.integral(1,2), 1./3*(pow(2,3)-1), epsilon, "g3");
    integral_value = fg.integral(0.1,10);
    limit = fg.integral_limit_b(0.1,integral_value);
    check_close(*limit, 10, epsilon, "g4");
    limit = fg.integral_limit_b(10,-integral_value);
    check_close(*limit, 0.1, 1e-8, "g4b");

    function<piecewise_power> fh{{0.1,0.5,1,7,8,12},{1,0.1,2,3,1,6}};
    integral_value = fh.integral(0.01,10);
    limit = fh.integral_limit_b(0.01,integral_value);
    check_close(*limit, 10, epsilon, "h1");
    limit = fh.integral_limit_b(10,-integral_value);
    check_close(*limit, 0.01, 1e-9, "h1b");
    function<piecewise_linear> fh_l{{0.1,0.5,1,7,8,12},{1,0.1,2,3,1,6}};
    check_small(*fh_l.integral_limit_b(10,-integral_value), 0.2, "h1b_l");
    
    pp_function fi{{1,2},{1,1}};
    pe_function fj{{1,2},{1,1}};
    pe_function fk;
    check_throw(fk.value(0));
    pe_function fl(2);
    check_close(fl.value(),2,1e-12,"fl");
    pp_function fm{{1, 2, 3, 4},{1,4,9,16}};
    fm.add_extrapolation_points(1);
    check_close(fm.value(0.01),1,1e-9);    
    fm.scale_x(2);
    check_close(fm.value(2),1,1e-12);
    fm.scale_y(2);
    check_close(fm.value(2),2,1e-12);
    pp_function fmb = significant_digits(fm, 10, 1);
    check_close(fmb.value(8),30,1e-12);
    check_throw(fmb.add_extrapolation_points(0));
    
    function<piecewise_linear> fn{{-2,-1,5},{-2,-1,5}};
    check_close(fn.value(-3),-3,1e-12,"fna");
    check_close(fn.derivative(0),1,1e-12,"fnb");
    check_close(fn.integral(-4,6),10,1e-12,"fnc");
    check_close(*fn.integral_limit_b(-4,10),6,1e-12,"fnd");
    check_close(*fn.integral_limit_b(6,-10),-4,1e-12,"fne");
    fn.add_extrapolation_points(0);
    check_small(fn.value(-5),1e-12,"fnf");
       
    function<piecewise_linear> fo{{-2,-1,5},{0,0,0}};
    check_small(fo.value(-3),1e-12,"foa");
    check_small(fo.derivative(0),1e-12,"fob");
    check_small(fo.integral(-4,6),1e-12,"foc");
    check(fo.integral_limit_b(-4,10).has_value()==false,"fod");

    function<piecewise_linear> fp{{-1,0,1},{0,1,2}};
    auto acc = fp.accumulation();
    check_close(acc.back(),2,1e-12,"fpa");
    auto ia = inverted_accumulation(fp);
    check_close(ia.value(2),1,1e-12,"fpb");
  } end_test_case()
}
