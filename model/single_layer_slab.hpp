#ifndef flick_single_layer_slab
#define flick_single_layer_slab

#include "../numeric/named_bounded_types.hpp"
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
    double relative_depth_{0};
    std::shared_ptr<transporter::ordinary_mc> omc_;
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
    double relative_radiance(const polar_angle& pa,
			     const azimuth_angle& aa,
			     const polar_angle& acceptance_angle
			     = polar_angle{10},
			     const unit_interval& relative_depth
			     = unit_interval{0}) {
      set_relative_depth(relative_depth());
      run();
      find_receivers();
      unit_vector direction{pa(),aa()};
      double dI_r = reflected_->projected_intensity(direction,
						    acceptance_angle());
      double dI_t = transmitted_->projected_intensity(direction,
						      acceptance_angle());
      return (dI_r + dI_t) / incident_->radiant_flux();
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
  private:
    void set_relative_depth(double rd) {
      relative_depth_ = rd;
      if (relative_depth_ < geometry_.small_step())
	relative_depth_ = geometry_.small_step()*2;
    }
    void find_receivers() {
      incident_ = &omc_->inward_receiver("surface");
      reflected_ = &omc_->outward_receiver("sheet");
      transmitted_ = &omc_->inward_receiver("sheet");
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
      unit_vector direction{constants::pi-theta_0_(),0}; 
      emitter.set_direction<unidirectional>(direction);
      build_geometry();
      double g = material_->asymmetry_factor();
      omc_ = std::make_shared<transporter::ordinary_mc>(geometry_);
      omc_->set(sampling_asymmetry_factor{g});
      omc_->transport_radiation(emitter,"geometry");
    }
  };
}
}

#endif
