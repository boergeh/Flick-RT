#include "pose.hpp"

namespace flick {
  using namespace constants;
  begin_test_case(pose_test) {
    pose p;
    p.rotate_about({1,0,0},pi);
    check_small(rms(p.z_direction(),{0,0,-1}));
    vector v = p.apply_pose_rotation_to({0,1,0});
    check_small(rms(v,{0,-1,0}));

    pose p1;
    pose p2;
    p1.move_to({3,0,0});
    p1.rotate_about({0,0,1},pi/2);
    p2.move_to({3,3,0});
    p2.rotate_about({1,0,0},pi/2);
    pose p3 = p2.as_observed_by(p1); 
    check_small(rms(p3.position(),{3,0,0}));
    check_small(rms(p3.x_direction(),{0,-1,0}));
    check_small(rms(p3.y_direction(),{0,0,1}));
    check_small(rms(p3.z_direction(),{-1,0,0}));

    pose observer({0,0,3},rotation_about_x(pi));
    pose placement{{0,0,2},no_rotation()};
    pose o = observer.as_observed_by(placement);
    check_small(rms(o.position(),{0,0,1}));
    check_small(rms(o.z_direction(),{0,0,-1}));

    unit_vector uv = {1,0,0};
    check_small(rms({-1,0,0},-uv));

    pose p4;
    unit_vector uv4 = {3*pi/2,0.1};
    p4.rotate_to(-uv4);
    check_small(dot(uv4,p4.z_direction())+1);
  } end_test_case()
}
