#ifndef flick_single_layer_slab
#define flick_single_layer_slab

#include "../material/named_bounded_types.hpp"
#include "../material/henyey_greenstein.hpp"
//#include "../geometry/volume.hpp"
//#include "../component/content.hpp"
//#include "../component/emitter.hpp"
#include "../transporter/ordinary_mc.hpp"

namespace flick {
namespace model {
  class single_layer_slab {
    thickness h_;
    incidence_angle theta_0_;
    bottom_albedo ba_;
    semi_infinite_box geometry_;
    number_of_packages np_;
    //precision precision_;
    std::shared_ptr<material::base> material_;
  public:
    single_layer_slab(const thickness& h) : h_{h} {
    }
    double hemispherical_reflectance() {
      // See Wikipedia reflectance for definition
      run();
      auto nav_ =  geometry::navigator<content>(geometry_);
      auto& v = nav_.find("layer");
      auto& inward = v().inward_receiver();
      auto& outward = v().outward_receiver();
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
      auto& incident = vl().inward_receiver();
      auto& vb = nav_.find("bottom");
      auto& transmitted = vb().inward_receiver();
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
    void set(const number_of_packages& np) {
      np_ = np;
    }
    void set(const incidence_angle& ia) {
      theta_0_ = ia;
    }
    // void set(precision& p) {
    //  recision_ = p;
    //}
    template <class Material, class... Args>
    void fill(Args... a) {
      material_ = std::make_shared<Material>(a...);
    }
    friend std::ostream& operator<<(std::ostream &os,
				    single_layer_slab& sls) {
      //os << "single-layer slab with " << sls.hg_
      //	 << ", bottom albedo " << sls.ba_() << ", and incidence angle "
      //	 << sls.theta_0_(); 
      return os;
    }
  private:
    void build_geometry() {
      semi_infinite_box layer;
      semi_infinite_box bottom;
      geometry_.name("geometry");
      layer.name("layer");
      bottom.name("bottom");
      layer().inward_receiver().activate();
      layer().outward_receiver().activate();
      layer().fill(material_);
      bottom().inward_receiver().activate();
      bottom().coat<coating::grey_lambert>(ba_(),1-ba_());
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
      transporter::ordinary_mc(geometry_,geometry_,emitter).run(1,0.5);
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
