#ifndef flick_single_layer_slab
#define flick_single_layer_slab

#include "named_bounded_types.hpp"
#include "../material/henyey_greenstein.hpp"
//#include "../geometry/volume.hpp"
//#include "../component/content.hpp"
//#include "../component/emitter.hpp"
#include "../transporter/ordinary_mc.hpp"

namespace flick {
namespace model {
  class single_layer_slab {
    thickness h_;
    material::henyey_greenstein hg_;
    incidence_angle theta_0_;
    bottom_albedo ba_;
    geometry::semi_infinite_box<content> geometry_;
    transporter::ordinary_mc transporter_;
  public:
    single_layer_slab(const thickness& h,
		      const material::henyey_greenstein& hg,
		      const incidence_angle& theta_0,
		      const transporter::ordinary_mc& tr)
      : h_{h}, hg_{hg}, theta_0_{theta_0}, transporter_{tr} {}
    double hemispherical_reflectance() {
      run(transporter_.np_());
      auto nav_ =  geometry::navigator<content>(geometry_);
      auto& v = nav_.find("layer");
      auto& inward = v.content().inward_receiver();
      auto& outward = v.content().outward_receiver();
      return outward.radiant_flux()/inward.radiant_flux();
    }
    void set(const bottom_albedo& ba) {
      ba_ = ba;
    }
    friend std::ostream& operator<<(std::ostream &os,
				    single_layer_slab& sls) {
      os << "single-layer slab with " << sls.hg_
	 << ", bottom albedo " << sls.ba_() << ", and incidence angle "
	 << sls.theta_0_(); 
      return os;
    }
  private:
    void build_geometry() {
      geometry::semi_infinite_box<content> layer;
      geometry::semi_infinite_box<content> bottom;
      geometry_.name("geometry");
      layer.name("layer");
      bottom.name("bottom");
      layer.content().inward_receiver().activate();
      layer.content().outward_receiver().activate();
      bottom.content().inward_receiver().activate();
      bottom.content().set_coating<coating::grey_lambert>(ba_(),1-ba_());
      geometry_.move_by({0,0,h_()+1});
      layer.move_by({0,0,h_()});
      layer.insert(bottom);
      geometry_.clear();
      geometry_.insert(layer);
    }
    void run(size_t n_packages) {
      emitter emitter{{0,0,h_()+0.5},stokes{1,0,0,0},n_packages};
      emitter.set_wavelength<monocromatic>(500e-9);
      auto direction = unit_vector{constants::pi-theta_0_(),0}; 
      emitter.set_direction<unidirectional>(direction);
      build_geometry();
      mc_basic(geometry_).run(emitter, geometry_);
    }
  };
}
}

#endif
