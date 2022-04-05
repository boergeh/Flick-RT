#include "pose.hpp"
namespace flick {
  begin_test_case(pose_test) {
    pose p;
    p.rotate_about({1,0,0},pi);
    check_small(rms(p.z_direction(),{0,0,-1}),1e-15);
    vector v = p.apply_pose_rotation_to({0,1,0});
    check_small(rms(v,{0,-1,0}),1e-15);

    pose p1;
    pose p2;
    p1.move_to({3,0,0});
    p1.rotate_about({0,0,1},pi/2);
    p2.move_to({3,3,0});
    p2.rotate_about({1,0,0},pi/2);
    pose p3 = p2.as_observed_by(p1); 
    check_small(rms(p3.position(),{3,0,0}),1e-15,"a");
    check_small(rms(p3.x_direction(),{0,-1,0}),1e-15,"b");
    check_small(rms(p3.y_direction(),{0,0,1}),1e-15,"c");
    check_small(rms(p3.z_direction(),{-1,0,0}),1e-15,"d");
    //show_position_and_directions(p3);

    pose observer({0,0,3},rotation_about_x(pi));
    pose placement{{0,0,2},no_rotation()};
    pose o = observer.as_observed_by(placement);
    check_small(rms(o.position(),{0,0,1}),1e-15);
    check_small(rms(o.z_direction(),{0,0,-1}),1e-15);
    //show_position_and_directions(op);

    //pose p_1 =  pose{{-1,0,0},no_rotation()};
    //pose p_2 =  pose{{1,0,0},rotation_about_x(pi)};
    //vector lp = local_point({0,0,1}).in(p_1).as_observed_by(p_2);
    //check_small(rms(lp,{-2,0,-1}),1e-15);
    //check_small(rms(2*lp,{-4,0,-2}),1e-15);

    unit_vector uv = {1,0,0};
    check_small(rms({-1,0,0},-uv),1e-15);

    pose p4;
    unit_vector uv4 = {3*pi/2,0.1};
    p4.rotate_to(-uv4);
    check_small(dot(uv4,p4.z_direction())+1,1e-15);

    //std::cout << std::endl;
    //pose p5{{1,0,0},rotation_to({1,0,0})};
    //std::cout << p5.position() << "  "<< p5.z_direction()<< std::endl;
    //pose p6 = p5.as_observed_by(pose{{-1,0,0}, rotation_to({1,0,0})});
    //std::cout << p6.position() << "  "<< p6.z_direction()<< std::endl;

    
  } end_test_case()
}
