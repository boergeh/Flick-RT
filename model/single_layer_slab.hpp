#ifndef flick_single_layer_slab
#define flick_single_layer_slab

#include "../material/named_bounded_types.hpp"
#include "../material/henyey_greenstein.hpp"
#include "../transporter/ordinary_mc.hpp"

namespace flick {
namespace model {
  class single_layer_slab
  // remote sensing reflectance method?
  // set precision instead of packages? 
  {
    thickness h_;
    incidence_angle theta_0_;
    bottom_albedo ba_;
    semi_infinite_box geometry_;
    number_of_packages np_;
    std::shared_ptr<material::base> material_;
    receiver* incident_;
    receiver* transmitted_;
    receiver* reflected_;
    double relative_depth_{0.5};
  public:
    single_layer_slab(const thickness& h) : h_{h} {}
    void set(const bottom_albedo& ba) {
      ba_ = ba;
    }
    void set(const number_of_packages& np) {
      np_ = np;
    }
    void set(const incidence_angle& ia) {
      theta_0_ = ia;
    }
    template <class Material, class... Args>
    void fill(Args... a) {
      material_ = std::make_shared<Material>(a...);
    }
    double hemispherical_reflectance(const unit_interval& relative_depth
				     = unit_interval{0})
    // See Wikipedia reflectance for definition
    {
      set_relative_depth(relative_depth());
      run();
      find_receivers();
      return reflected_->radiant_flux()/incident_->radiant_flux();
    }
    double directional_reflectance(const polar_angle& pa,
				   const azimuth_angle& aa,
				   const unit_interval& relative_depth
				     = unit_interval{0})
    // See Wikipedia reflectance for definition

    // should search for optimal angle. Gaussian aceptance?
    // should guess direction directly towards detector?
    // should sample with high prob. at detector depth?
    {
      set_relative_depth(relative_depth());
      run();
      find_receivers();
      unit_vector direction{pa(),aa()};
      return reflected_->radiant_intensity(direction,20/180*constants::pi)/incident_->radiant_flux();
    }
    double hemispherical_transmittance(const unit_interval& relative_depth
				       = unit_interval{1})
    // See Wikipedia transmittance for definition
    {
      set_relative_depth(relative_depth());
      run();
      find_receivers();
      return transmitted_->radiant_flux()/incident_->radiant_flux();
    }
    double directional_transmittance(const polar_angle& pa,
				     const azimuth_angle& aa,
				     const unit_interval& relative_depth
				     = unit_interval{0})
    // See Wikipedia transmittance for definition
    {
      set_relative_depth(relative_depth());
      run();
      find_receivers();
      unit_vector direction{pa(),aa()};
      return transmitted_->radiant_intensity(direction,20/180*constants::pi)/incident_->radiant_flux();
    }
  private:
    void set_relative_depth(double rd) {
      relative_depth_ = rd;
      if (relative_depth_ < geometry_.small_step())
	relative_depth_ = geometry_.small_step()*2;
    }
    void find_receivers() {
      auto nav =  geometry::navigator<content>(geometry_);
      auto& v_in = nav.find("surface");
      incident_ = &v_in().inward_receiver();
      auto& v_sheet = nav.find("sheet");
      reflected_ = &v_sheet().outward_receiver();
      transmitted_ = &v_sheet().inward_receiver();
    }
    void build_geometry() {
      semi_infinite_box surface;
      semi_infinite_box sheet;
      semi_infinite_box bottom;
      geometry_.name("geometry");
      surface.name("surface");
      sheet.name("sheet");
      bottom.name("bottom");
      surface().inward_receiver().activate();
      surface().fill(material_);
      sheet().inward_receiver().activate();
      sheet().outward_receiver().activate();
      sheet().fill(material_);
      bottom().coat<coating::grey_lambert>(ba_(),1-ba_());
      geometry_.move_by({0,0,h_()+1});
      surface.move_by({0,0,h_()});
      sheet.move_by({0,0,h_()*(1-relative_depth_)});
      sheet.insert(bottom);
      surface.insert(sheet);
      geometry_.clear();
      geometry_.insert(surface);
    }
    void run() {
      emitter emitter{{0,0,h_()+0.5},stokes{1,0,0,0},np_()};
      emitter.set_wavelength<monocromatic>(500e-9);
      auto direction = unit_vector{constants::pi-theta_0_(),0}; 
      emitter.set_direction<unidirectional>(direction);
      build_geometry();
      double g = material_->sampling_asymmetry_factor();
      transporter::ordinary_mc(geometry_,geometry_,emitter).run(1,g);
    }
  };
}
}

#endif
