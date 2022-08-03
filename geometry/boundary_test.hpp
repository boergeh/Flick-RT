#include "boundary.hpp"

namespace flick {
  begin_test_case(boundary_test) {
    using namespace constants;
    using namespace geometry; 
    using namespace geometry::surface;
    pose observer({0,0,3},rotation_about_x(pi));
    auto b = boundary().add(make_surface<plane>(),pose{{0,0,2}, no_rotation()});
    vector p = (*b.intersection(observer)).position();
    check_small(rms(p,{0,0,2}),1e-15,"a");
    pose observer2({0,0,3},rotation_about_x(3*pi/2));
    vector surf_normal = (*b.intersection(observer2)).z_direction();
    check_small(rms(surf_normal,{0,0,1}),1e-15,"b");
    
    b.add(make_surface<sphere>(1));
    vector pos = (*b.intersection(observer)).position();
    vector n = (*b.intersection(observer)).z_direction();
    check_small(rms(pos,{0,0,1}),1e-15,"c");
    check_small(rms(n,{0,0,1}),1e-15,"d");
    p = (*b.intersection(observer.move_to({0,0,0}))).position();
    check_small(rms(p,{0,0,-1}),1e-15,"e");
    p = (*b.intersection(observer.rotate_about(observer.x_direction(),pi))).position();
    check_small(rms(p,{0,0,1}),1e-15,"f");
    observer = pose{{0,0,0},rotation_about_y(pi/4)};
    p = (*b.intersection(observer)).position();
    check_small(rms(p,{1/sqrt(2),0,1/sqrt(2)}),1e-15,"g");
    observer.move_to({3,0,0}).rotate_to({1,0,0});
    
    b.add(make_surface<sphere>(0.5));
    p = (*b.intersection({{0,0,-10}, no_rotation()})).position();
    check_small(rms(p,{0,0,-0.5}),1e-15,"h");
    
    b.add(make_surface<sphere>(0.25),default_pose,inside_out);
    p = (*b.intersection({{0,0,-10}, no_rotation()})).position();
    check_small(rms(p,{0,0,-0.5}),1e-15,"i");
    
    vector center{0,1,2};
    boundary b2;
    b2.add(make_surface<sphere>(1));
    b2.move_by(center);
    auto ui2 = uniform_intersections{b2,3000,limits{0.1,12}};
    check_close(ui2.boundary_area(),4*pi, 9, "j");
    check_close(ui2.enclosed_volume(),4./3*pi, 9, "k");
    check_small(rms(ui2.center_of_gravity(),center),0.2,"l");
    
    boundary b3;
    double r1 = 2;
    b3.add(make_surface<sphere>(r1));
    b3.add(make_surface<plane>(),{{0,0,0},rotation_to({0,0,-1})});
    pose observer3{{3,0,-1},no_rotation()};
    check(b3.intersection(observer3).has_value()==false,"m");
    check(b3.intersection(observer3.move_to({0,0,-9})).has_value()==true,"n");
    check(b3.intersection(observer3.move_to({0,0,-0.1})).has_value()==true,"o");
    auto ui3 = uniform_intersections{b3,2000,limits{3,3.2}};
    check_close(ui3.enclosed_volume(),4./3*pi*pow(r1,3)/2, 9, "p");
    double a = 4*pi*pow(r1,2)/2+pi*pow(r1,2);
    check_close(ui3.boundary_area(),a, 9, "q");
    
    boundary b4;
    double d = 15;
    auto pl4 = std::make_shared<surface::plane>();
    vector center4 = {30,30,30};
    b4.add(pl4,{{0,0,d}, rotation_to({0,0,1})});
    b4.add(pl4,{{0,0,-d}, rotation_to({0,0,-1})});
    b4.add(pl4,{{0,d,0}, rotation_to({0,1,0})});
    b4.add(pl4,{{0,-d,0}, rotation_to({0,-1,0})});
    b4.add(pl4,{{d,0,0}, rotation_to({1,0,0})});
    b4.add(pl4,{{-d,0,0}, rotation_to({-1,0,0})});
    b4.move_by(center4);
    b4.rotate_by(rotation_about_x(pi/4),b4.placement().position());
    auto ui4 = uniform_intersections{b4,10000,limits{1.0,40}};
    check_close(ui4.boundary_area(),pow(2*d,2)*6,9,"r");
    check_close(ui4.enclosed_volume(),pow(2*d,3),9,"s");
    check_small(rms(ui4.center_of_gravity(),center4),2,"t");

    boundary b5;
    auto sp5 = make_surface<sphere>(1);
    auto pl5 = make_surface<plane>();
    b5.add(sp5);    
    b5.add(pl5);
    b5.add(pl5,{{0,0,0}, rotation_about_y(pi/2)});
    b5.move_by({0,0,2});
    b5.rotate_by(rotation_about_y(pi/9));
    auto ui5 = uniform_intersections{b5,10000,limits{1.1,20}};  
    check_close(ui5.boundary_area(),4*pi/4+pi,9,"u");
    check_close(ui5.enclosed_volume(),4./3*pi/4,9,"v");
   
    boundary b6;
    auto sp6 = std::make_shared<surface::sphere>(1);
    b6.add(sp6);
    b6.add(make_surface<sphere>(0.5),default_pose,inside_out);
    auto ui6 = uniform_intersections{b6,2000,limits{1.1,20}};  
    check_close(ui6.boundary_area(),4*pi*(1+pow(0.5,2)),9,"w");
    check_close(ui6.enclosed_volume(),4./3*pi*(1-pow(0.5,3)),9,"x");

    //write<uniform_intersections>(ui6,"boundary_view.txt");

    cubical_boundary cb{1};
    cb.move_by({1,0,0});
    cb.rotate_by(rotation_about_z(pi/2),{-1,0,0});
    check_small(rms(cb.placement().position(),{-1,2,0}),1e-9);
    
  } end_test_case()
}
