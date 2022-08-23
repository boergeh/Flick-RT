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
    number_of_packages np_;
  public:
    single_layer_slab(const thickness& h,
		      const material::henyey_greenstein& hg,
		      const incidence_angle& theta_0,
		      const number_of_packages& np)
      : h_{h}, hg_{hg}, theta_0_{theta_0}, np_{np} {}
    double hemispherical_reflectance() {
      // See Wikipedia reflectance for definition
      run();
      auto nav_ =  geometry::navigator<content>(geometry_);
      auto& v = nav_.find("layer");
      auto& inward = v.content().inward_receiver();
      auto& outward = v.content().outward_receiver();
      return outward.radiant_flux()/inward.radiant_flux();
    }
    double directional_reflectance(double polar_angle,
				   double azimuth_angle)
    // See Wikipedia reflectance for definition
    {
      //tbi
      // should search for optimal angle. Gaussian aceptance?
      // should gess direction directly towards detector
      return 0;
    }
    double hemispherical_transmittance()
    // See Wikipedia transmittance for definition
    {
      run();
      auto nav_ =  geometry::navigator<content>(geometry_);
      auto& vl = nav_.find("layer");
      auto& incident = vl.content().inward_receiver();
      auto& vb = nav_.find("bottom");
      auto& transmitted = vb.content().inward_receiver();
      return transmitted.radiant_flux()/incident.radiant_flux();
    }
    double directional_transmittance(double polar_angle,
				     double azimuth_angle)
    // See Wikipedia transmittance for definition
    {
      //tbi
      return 0;
    }
    void set(const bottom_albedo& ba) {
      ba_ = ba;
    }
    void set(number_of_packages& np) {
      np_=np;
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
    void run() {
      emitter emitter{{0,0,h_()+0.5},stokes{1,0,0,0},np_()};
      emitter.set_wavelength<monocromatic>(500e-9);
      auto direction = unit_vector{constants::pi-theta_0_(),0}; 
      emitter.set_direction<unidirectional>(direction);
      build_geometry();
      transporter::ordinary_mc(geometry_,geometry_,emitter).run();
    }
  };
  /*
   double remote_sensing_reflectance(const basic_slab& bs) {
	//tbi
	return 0;
      }
      class henyey_greenstein_slab : public basic_slab {
	double asymmetr_factor_{0};
      public:	
      };
      class user_defined_slab : public basic_slab {
	
      public:
	void tabulated_phase_function(const std::string& filename) {
	}
      };
  */
}
}

#endif
