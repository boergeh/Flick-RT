#include "volume.hpp"
namespace flick {
namespace geometry {
  begin_test_case(volume_test) {
    class content
    // User defined volume content. May include material, coating,
    // etc.
    {      
    public:
      double mass{9};
    };
    cube<content> c{1};       
    auto ui = c.get_uniform_intersections(10000);
    //write<uniform_intersections>(ui,"volume_view.txt",4);
    check_close(ui.boundary_area(),6,9);
    check_close(ui.enclosed_volume(),1,9);
    check_small(rms(ui.center_of_gravity(),{0,0,0}),0.1);

    cube<content> c2{1};
    check_small(c2().mass-9,1e-9);
  
    sphere<content> small_sphere{2};
    small_sphere.insert(c2);
     
    sphere<content> large_sphere{3};
    large_sphere.insert(small_sphere);
    navigator<content> nav(large_sphere);
    nav.go_inward();
    nav.go_inward();
    nav.go_outward();
    nav.go_outward();
    nav.go_to("cube");
    nav.go_outward();
    nav.go_outward(); 
    pose observer{{0,0,-10},{0,0}};
    check_throw(nav.next_intersection(observer));
    observer.move_to({0,0,-2.99});
    vector pos = (*nav.next_intersection(observer)).position();
    check_small(rms(pos,{0,0,-2}),1e-9);
    check(nav.next_volume(observer).name()=="sphere");
    nav.go_inward();
    observer.move_to({0,0,-1.99});
    pos = (*nav.next_intersection(observer)).position();
    check_small(rms(pos,{0,0,-0.5}),1e-9);
    check(nav.next_volume(observer).name()=="cube");
    nav.go_inward();
    check_throw(nav.go_inward());
    nav.go_outward();
    nav.go_outward();
    check_throw(nav.go_outward());
    
    sphere sph_copy = large_sphere;
    vector new_pos = {1,1,1};
    sph_copy.move_by(new_pos);
    navigator nav2(sph_copy);
    nav2.go_inward();
    nav2.go_inward();
    vector sph_pos = nav2.go_outward().placement().position();
    check_small(rms(sph_pos,new_pos),1e-9);

    cube<content> c3{0.1};
    c3.move_by({1,0,0});
    c3.rotate_by(rotation_about_z(constants::pi/2),{-1,0,0});
    check_small(rms(c3.placement().position(),{-1,2,0}),1e-9);
    
  } end_test_case()
}
}
